#ifndef _EXCHANGEBUFFER_HPP_
#define _EXCHANGEBUFFER_HPP_

#include "Buffer.hpp"

#include "Logging.hpp"

#include <vector>

class ExchangeBuffer {
private:
  std::vector<int> numOfMessages;

  std::vector<int> destinations;

  std::vector<size_t> sendSizes;

  int numberOfSendMessages;

  int numberOfReceiveMessages;

  bool init;

  int commSplit;

  std::vector<Buffer> sendBuffers;

  std::vector<Buffer> receiveBuffers;

  void CommunicateSizes();
public:
  bool IsInitialized() const;

  void CommunicateSize();

  size_t ReceiveSize(int q) const;

  void SendSize(size_t m, int q);

  int SendMessages() const;

  int RecvMessages() const;

  int Messages(int q) const;

  int MessageDest(int q) const;

  size_t MessageSize(int q) const;

  friend std::ostream &operator<<(std::ostream &os, const ExchangeBuffer &exBuffer);

  ExchangeBuffer(int commSplit = 0);

  ~ExchangeBuffer();

  void Destruct();

  // Every proc calling this function sends TO proc q
  Buffer &Send(int q);

  // Only proc which has been addressed in Send receives message
  Buffer &Receive(int q);

  void Communicate();

  void ClearBuffers();

  ExchangeBuffer &Rewind();

  int CommSplit() const;
};

#endif // of #ifndef _EXCHANGEBUFFER_HPP_
