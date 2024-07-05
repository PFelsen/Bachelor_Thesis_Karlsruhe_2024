#!/bin/bash

cd ..

WORKSPACE=$(pwd)

echo "Use Mpp Data Git: $USE_MPP_DATA"
echo "Use LSDF: $USE_LSDF"

if [ $USE_LSDF == "true" ]; then
  LSDF_DIR=$LSDFPROJECTS/mppci/$CI_PROJECT_NAME/$CI_JOB_NAME/$(date +"%Y-%m-%d_%T")
  echo "Store in LSDF with path:"
  echo $LSDF_DIR
  mkdir -p $LSDF_DIR
  cd $LSDF_DIR
  mkdir vtu
  mkdir py
  mkdir vtk
  mkdir log
else
  echo "Don't store in LSDF"
fi

if [ -z "$(ls -A $WORKSPACE/data/vtu/)" ]; then
  echo "No vtu files found"
else
  cp $WORKSPACE/data/vtu/*.vtu $CI_PROJECT_DIR
  if [ $USE_LSDF == "true" ]; then
    cp $WORKSPACE/data/vtu/*.vtu $LSDF_DIR/vtu/
  fi
  rm $WORKSPACE/data/vtu/*.vtu
fi

if [ -z "$(ls -A $WORKSPACE/json/)" ]; then
  echo "No json files found"
else
  cp $WORKSPACE/json/*.json $CI_PROJECT_DIR
  if [ $USE_LSDF == "true" ]; then
    cp $WORKSPACE/json/*.json $LSDF_DIR/json/
  fi
  rm $WORKSPACE/json/*.json
fi


if [ -z "$(ls -A $WORKSPACE/data/py/)" ]; then
  echo "No png files found"
else
  cp $WORKSPACE/data/py/*.png $CI_PROJECT_DIR
  if [ $USE_LSDF == "true" ]; then
    cp $WORKSPACE/data/py/*.png $LSDF_DIR/png/
  fi
  rm $WORKSPACE/data/py/*.png
fi


if [ -z "$(ls -A $WORKSPACE/log/)" ]; then
  echo "No log files found"
else
  cp $WORKSPACE/log/*.log $CI_PROJECT_DIR
  if [ $USE_LSDF == "true" ]; then
    cp $WORKSPACE/log/*.log $LSDF_DIR/log/
  fi
  rm $WORKSPACE/log/*.log
fi


if [ -z "$(ls -A $WORKSPACE/data/vtk/)" ]; then
  echo "No vtk files found"
else
  cp $WORKSPACE/data/vtk/*.vtk $CI_PROJECT_DIR
  if [ $USE_LSDF == "true" ]; then
    cp $WORKSPACE/data/vtk/*.vtk $LSDF_DIR/vtk/
  fi
  rm $WORKSPACE/data/vtk/*.vtk
fi


if [ $USE_MPP_DATA == "true" ]; then
  git config --global user.name "${GITLAB_USER_NAME}"
  git config --global user.email "${GITLAB_USER_EMAIL}"
  rm -rf mpp-data
  git clone https://gitlab-ci-token:${MPP_DATA_ACCESS}@$CI_SERVER_HOST/kit/mpp/mpp-data.git
  cd mpp-data
  git status
  git pull origin main
  mkdir -p $CI_PROJECT_NAME/$CI_JOB_NAME/json
  mkdir -p $CI_PROJECT_NAME/$CI_JOB_NAME/vtu
  mkdir -p $CI_PROJECT_NAME/$CI_JOB_NAME/vtk
  mkdir -p $CI_PROJECT_NAME/$CI_JOB_NAME/py
  mkdir -p $CI_PROJECT_NAME/$CI_JOB_NAME/log
  cp $CI_PROJECT_DIR/*.vtk $CI_PROJECT_NAME/$CI_JOB_NAME/vtk/
  cp $CI_PROJECT_DIR/*.log $CI_PROJECT_NAME/$CI_JOB_NAME/log/
  cp $CI_PROJECT_DIR/*.png $CI_PROJECT_NAME/$CI_JOB_NAME/py/
  cp $CI_PROJECT_DIR/*.vtu $CI_PROJECT_NAME/$CI_JOB_NAME/vtu/
  cp $CI_PROJECT_DIR/*.json $CI_PROJECT_NAME/$CI_JOB_NAME/json/
  git add -f --all
  git commit -m "Updated data for $BENCHMARK at $(date)"
  git push https://gitlab-ci-token:${MPP_DATA_ACCESS}@$CI_SERVER_HOST/kit/mpp/mpp-data.git HEAD:main
fi

