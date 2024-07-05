#include "TestEnvironment.hpp"

#include "Buffer.hpp"

class TestBuffer : public Test {
protected:
  Buffer buffer;

  void TearDown() override { buffer.Destruct(); }
};

TEST_F(TestBuffer, TestBufferChar) {
  //    const char send = 'A';
  //    char recv = 'B';
  //    buffer.fill(send, 1);
  //    buffer.read(recv, 1);
  //    EXPECT_TRUE(send == recv);
}

TEST_F(TestBuffer, TestBufferChar2) {
  const char sendChar = 'A';
  const char *send = &sendChar;
  char recvChar = 'B';
  char *recv = &recvChar;
  buffer << send;

  std::vector<int> vec(2);
  buffer >> recv;
  EXPECT_TRUE(*send == *recv);
}

TEST_F(TestBuffer, TestBufferString) {
  //    const std::string send = "Test";
  //    std::string recv;
  //    buffer.fill(send, sizeof(char));
  //    buffer.read(recv, sizeof(char));
  //    EXPECT_TRUE(send == recv);
}

TEST_F(TestBuffer, TestBufferString2) {
  //    const std::string send = "Test";
  //    std::string recv;
  //    buffer.resize(send.length());
  //    for (int i = 0; i < send.length(); i++)
  //    buffer << send;
  //    buffer >> recv;
  //    EXPECT_TRUE(send == recv);
}

int main(int argc, char **argv) { return MppTest(MppTestBuilder(argc, argv)).RUN_ALL_MPP_TESTS(); }
