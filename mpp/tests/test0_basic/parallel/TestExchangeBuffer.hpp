#ifndef TESTEXCHANGEBUFFER_HPP
#define TESTEXCHANGEBUFFER_HPP

#include <array>
#include <list>
#include <set>
#include <string>
#include <tuple>
#include <vector>
#include "TestEnvironment.hpp"

// Static Tests for concepts
static_assert(MapIterableLike<std::unordered_map<int, int>>);
static_assert(MapIterableLike<std::map<int, int>>);
static_assert(MapIterableLike<std::multimap<int, int>>);

static_assert(SizableIterableLike<std::vector<int>>);
static_assert(SizableIterableLike<std::set<int>>);
static_assert(SizableIterableLike<std::list<int>>);
static_assert(SizableIterableLike<std::multiset<int, int>>);
static_assert(SizableIterableLike<std::string>);

static_assert(ReservableContainer<std::vector<int>>);
static_assert(ReservableContainer<std::unordered_map<int, int>>);
static_assert(ReservableContainer<std::string>);

static_assert(TupleLike<std::tuple<int, int, int, int>>);
static_assert(TupleLike<std::pair<int, int>>);
static_assert(TupleLike<std::tuple<int>>);
static_assert(TupleLike<std::array<int, 4>>);

// End static asserts

/*
 * Improvements:
 *  - it would be much nicer if every proc check is its own test case
 *  - also, tests should focus more on individual functions instead of
 *      overall class functionality
 */

class TestExchangeBuffer : public Test {
protected:
  int commSplit = 0;

  ExchangeBuffer *exBuffer;

  void SetUp() override {
    PPM->Barrier(commSplit);
    exBuffer = new ExchangeBuffer(0);
  }

  void TearDown() override {
    PPM->Barrier(commSplit);
    delete exBuffer;
  }

  void checkSize(int q, size_t sizeOnProc) {
    if (PPM->Proc(commSplit) == q) {
      for (int p = 0; p < PPM->Size(commSplit); p++)
        EXPECT_EQ(exBuffer->ReceiveSize(p), sizeOnProc);
    } else {
      EXPECT_EQ(exBuffer->ReceiveSize(q), 0);
    }
  }

  template<typename T>
  void checkMessage(int q, T send, T recv) {
    if (PPM->Proc(commSplit) == q) {
      for (short p = 0; p < PPM->Size(commSplit); ++p) {
        while (exBuffer->Receive(p).size() < exBuffer->ReceiveSize(p)) {
          exBuffer->Receive(p) >> recv;
        }
        EXPECT_EQ(send, recv);
      }
    } else {
      EXPECT_NE(send, recv);
    }
  }

  void checkStrMessage(int q, const std::string &send) {
    std::string recv = "";
    if (PPM->Proc(commSplit) == q) {
      for (short p = 0; p < PPM->Size(commSplit); ++p) {
        while (exBuffer->Receive(p).size() < exBuffer->ReceiveSize(p)) {
          char character;
          exBuffer->Receive(p) >> character;
          recv += character;
        }
        EXPECT_EQ(send, recv);
        recv = "";
      }
    } else {
      EXPECT_NE(send, recv);
    }
  }

  void checkCStrMessage(int q, const char *send) {
    std::string recv = "";
    if (PPM->Proc(commSplit) == q) {
      for (short p = 0; p < PPM->Size(commSplit); ++p) {
        while (exBuffer->Receive(p).size() < exBuffer->ReceiveSize(p)) {
          char character;
          exBuffer->Receive(p) >> character;
          recv += character;
        }
        EXPECT_EQ(std::string(send), recv);
        recv = "";
      }
    } else {
      EXPECT_NE(std::string(send), recv);
    }
  }

  void checkPointMessage(int q, Point send, Point recv) {
    if (PPM->Proc(commSplit) == q) {
      for (short p = 0; p < PPM->Size(commSplit); ++p) {
        while (exBuffer->Receive(p).size() < exBuffer->ReceiveSize(p)) {
          exBuffer->Receive(p) >> recv;
        }
        EXPECT_POINT_EQ(send, recv);
      }
    } else {
      EXPECT_POINT_NE(send, recv);
    }
  }

  virtual void checkIsInitialized() { EXPECT_FALSE(exBuffer->IsInitialized()); }

  virtual void checkIsInitialized2() {
    exBuffer->CommunicateSize();
    EXPECT_TRUE(exBuffer->IsInitialized());
  }

  virtual void checkReceiveSize() { EXPECT_EQ(exBuffer->ReceiveSize(PPM->Proc(commSplit)), 0); }

  virtual void checkSendSize() {
    exBuffer->SendSize(1, PPM->Proc(commSplit));
    exBuffer->CommunicateSize();
    EXPECT_EQ(exBuffer->ReceiveSize(PPM->Proc(commSplit)), 1);
  }

  virtual void checkShort() {
    short send = 0;
    short recv = 1;
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << send;
      exBuffer->Communicate();
      checkSize(p, sizeof(short));
      checkMessage(p, send, recv);
      TearDown();
      SetUp();
    }
  }

  virtual void checkChar() {
    char send = 'A';
    char recv = 'B';
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << send;
      exBuffer->Communicate();
      checkSize(p, sizeof(char));
      checkMessage(p, send, recv);
      TearDown();
      SetUp();
    }
  }

  virtual void checkInt() {
    int send = 0;
    int recv = 1;
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << send;
      exBuffer->Communicate();
      checkSize(p, sizeof(int));
      checkMessage(p, send, recv);
      TearDown();
      SetUp();
    }
  }

  virtual void checkDouble() {
    double send = 0.0;
    double recv = 1.0;
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << send;
      exBuffer->Communicate();
      checkSize(p, sizeof(double));
      checkMessage(p, send, recv);
      TearDown();
      SetUp();
    }
  }

  virtual void checkCStr() {
    std::string test = "Test";
    const char *send = test.c_str();
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      for (int i = 0; i < test.length(); ++i)
        exBuffer->Send(p) << send[i];
      exBuffer->Communicate();
      checkSize(p, 4);
      checkCStrMessage(p, send);
      TearDown();
      SetUp();
    }
  }

  virtual void checkStr() {
    std::string send = "Test";
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      for (int i = 0; i < send.length(); ++i)
        exBuffer->Send(p) << send[i];
      exBuffer->Communicate();
      checkSize(p, 4);
      checkStrMessage(p, send);
      TearDown();
      SetUp();
    }
  }

  virtual void checkStr2() {
    std::string send = "A";
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      for (int i = 0; i < send.length(); ++i)
        exBuffer->Send(p) << send[i];
      exBuffer->Communicate();
      checkSize(p, 1);
      checkStrMessage(p, send);
      TearDown();
      SetUp();
    }
  }

  virtual void checkPoint() {
    Point send(0.0, 0, 0, 0);
    Point recv(1.0, 0.0, 0.0);
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << send;
      exBuffer->Communicate();
      checkSize(p, sizeof(send));
      checkPointMessage(p, send, recv);
      TearDown();
      SetUp();
    }
  }

  void checkUnorderedMessage(int q, std::unordered_map<Point, double> send,
                             std::unordered_map<Point, double> receive) {
    if (PPM->Proc(commSplit) == q) {
      for (short p = 0; p < PPM->Size(commSplit); ++p) {
        while (exBuffer->Receive(p).size() < exBuffer->ReceiveSize(p)) {
          exBuffer->Receive(p) >> receive;
        }
        EXPECT_EQ(send.size(), receive.size());
        for (const auto &[key, value] : send) {
          const auto findIterator = receive.find(key);
          EXPECT_NE(findIterator, std::end(receive));
          EXPECT_EQ(findIterator->second, value);
        }
      }
    } else {
      EXPECT_TRUE(receive.empty());
    }
  }

  virtual void checkUnorderedMap() {
    std::unordered_map<Point, double> sendMap;
    sendMap[Point(1.0, 1.0, 0.0)] = 93248720934.04234234;
    sendMap[Point(1.0, 13424.0, 0.0)] = 435345345345345;
    sendMap[Point(1.0, 1.0, 234234.0)] = 34534535345353453;
    sendMap[Point(1234234.0, 1.0, 0.0)] = 3453534534534534534;
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << sendMap;
      exBuffer->Communicate();
      checkSize(p, sendMap.size() * (sizeof(Point) + sizeof(double)) + sizeof(sendMap.size()));
      std::unordered_map<Point, double> receiveMap;
      checkUnorderedMessage(p, sendMap, receiveMap);
      TearDown();
      SetUp();
    }
  }

  void checkTupleMessage(int q, std::tuple<Point, double, double> &send,
                         std::tuple<Point, double, double> &receive) {
    if (PPM->Proc(commSplit) == q) {
      for (short p = 0; p < PPM->Size(commSplit); ++p) {
        while (exBuffer->Receive(p).size() < exBuffer->ReceiveSize(p)) {
          exBuffer->Receive(p) >> receive;
        }
        const auto &[sendPoint, sendD1, sendD2] = send;
        const auto &[point, d1, d2] = receive;
        EXPECT_POINT_EQ(sendPoint, point);
        EXPECT_EQ(sendD1, d1);
        EXPECT_EQ(sendD2, d2);
      }
    } else {
      const auto &[sendPoint, sendD1, sendD2] = send;
      const auto &[point, d1, d2] = receive;
      EXPECT_POINT_NE(sendPoint, point);
      EXPECT_NE(sendD2, d2);
      EXPECT_NE(sendD2, d2);
    }
  }

  virtual void checkTuple() {
    std::tuple<Point, double, double> sendMap{Point(1.0, 1.0, 0.0), 93248720934.04234234,
                                              93248720934.04234234};
    for (int p = 0; p < PPM->Size(commSplit); p++) {
      exBuffer->Send(p) << sendMap;
      exBuffer->Communicate();
      checkSize(p, sizeof(Point) + sizeof(double) * 2);
      std::tuple<Point, double, double> receiveMap;
      checkTupleMessage(p, sendMap, receiveMap);
      TearDown();
      SetUp();
    }
  }
};

class TestExchangeBufferWithSplit : public TestExchangeBuffer {
protected:
  void SetUp() override {
    commSplit = 1;
    PPM->SplitCommunicators();
    PPM->Barrier(commSplit);
    exBuffer = new ExchangeBuffer(commSplit);
  }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(commSplit);
    delete exBuffer;
  }
};

class TestExchangeBufferWithDoubleSplit : public TestExchangeBuffer {
protected:
  void SetUp() override {
    commSplit = 2;
    PPM->SplitCommunicators();
    PPM->SplitCommunicators();
    PPM->Barrier(commSplit);
    exBuffer = new ExchangeBuffer(commSplit);
  }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(commSplit);
    delete exBuffer;
  }
};

class TestExchangeBufferWithFullSplit : public TestExchangeBuffer {
protected:
  void SetUp() override {
    PPM->FullSplit();
    commSplit = PPM->MaxCommSplit();
    PPM->Barrier(commSplit);
    exBuffer = new ExchangeBuffer(commSplit);
  }

  void TearDown() override {
    PPM->ClearCommunicators();
    PPM->Barrier(commSplit);
    delete exBuffer;
  }
};

#define EXCHANGEBUFFER_TESTS(ExchangeBufferTestClass)                                              \
  TEST_F(ExchangeBufferTestClass, TestIsInitialized) {                                             \
    checkIsInitialized();                                                                          \
    checkIsInitialized2();                                                                         \
  }                                                                                                \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestReceiveSize) { checkReceiveSize(); }                         \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestSendSize) { checkSendSize(); }                               \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestShort) { checkShort(); }                                     \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestChar) { checkChar(); }                                       \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestInt) { checkInt(); }                                         \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestDouble) { checkDouble(); }                                   \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestCStr) { checkCStr(); }                                       \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestStr) { checkStr(); }                                         \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestStr2) { checkStr2(); }                                       \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestPoint) { checkPoint(); }                                     \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestUnorderedMap) { checkUnorderedMap(); }                       \
                                                                                                   \
  TEST_F(ExchangeBufferTestClass, TestTuple) { checkTuple(); }

#endif // TESTEXCHANGEBUFFER_HPP
