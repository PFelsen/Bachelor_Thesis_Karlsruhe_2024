#!/bin/bash
#Adds a prepare-commit-msg hook to mpp
#Purpose: Add branch name to commit message automatically
#Note that master branch is excluded

echo "-- Setting prepare commit message hook to mpp"

DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"

git=$DIR"/../.git"
if [ -f "$git" ]; then
  while read line; do
    id=${line:0:7}
    if [ "$id" = "gitdir:" ]; then
      gitpath=$DIR"/../"${line##*:}
    fi
  done < $git
else
  gitpath=$DIR"/../.git"
fi

gitpath="$(echo -e "${gitpath}" | tr -d '[:space:]')"

hookpath=$gitpath"/hooks/"
filenameCommitMsg="prepare-commit-msg"

echo "   file: "$hookpath$filenameCommitMsg
cd $hookpath

addToFile=true
if [ -f "$filenameCommitMsg" ]; then
  while read line; do
    if [ "$line" = "#Add branch name to commit message" ]; then
      addToFile=false
    fi
  done < $filenameCommitMsg
else
  touch $filenameCommitMsg
  echo "#!/bin/bash" >> $filenameCommitMsg
fi

if [ "$addToFile" = true ] ; then
  echo "" >> $filenameCommitMsg
  echo "#Add branch name to commit message" >> $filenameCommitMsg
  echo "if [ -z \"\$BRANCHES_TO_SKIP\" ]; then" >> $filenameCommitMsg
  echo "  BRANCHES_TO_SKIP=(master)" >> $filenameCommitMsg
  echo "fi" >> $filenameCommitMsg
  echo "" >> $filenameCommitMsg
  echo "BRANCH_NAME=\$(git symbolic-ref --short HEAD)" >> $filenameCommitMsg
  echo "BRANCH_NAME=\"\${BRANCH_NAME##*/}\"" >> $filenameCommitMsg
  echo "" >> $filenameCommitMsg
  echo "BRANCH_EXCLUDED=\$(printf \"%s\n\" \"\${BRANCHES_TO_SKIP[@]}\" | grep -c \"^$BRANCH_NAME\$\")" >> $filenameCommitMsg
  echo "BRANCH_IN_COMMIT=\$(grep -c \"\[\$BRANCH_NAME\]\" \$1)" >> $filenameCommitMsg
  echo "" >> $filenameCommitMsg
  echo "if [ -n \"\$BRANCH_NAME\" ] && ! [[ \$BRANCH_EXCLUDED -eq 1 ]] && ! [[ \$BRANCH_IN_COMMIT -ge 1 ]]; then" >> $filenameCommitMsg
  echo "  sed -i.bak -e \"1s/^/[\$BRANCH_NAME] /\" \$1" >> $filenameCommitMsg
  echo "fi"  >> $filenameCommitMsg

  chmod +x $filenameCommitMsg
fi




