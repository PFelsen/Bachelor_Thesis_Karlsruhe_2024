#include "ExchangeBuffer.hpp"
#include "Parallel.hpp"

ExchangeBuffer::ExchangeBuffer(int commSplit) :
    numOfMessages(PPM->Size(commSplit) + 1), sendSizes(PPM->Size(commSplit), 0),
    numberOfSendMessages(0), numberOfReceiveMessages(0), init(false), commSplit(commSplit),
    sendBuffers(PPM->Size(commSplit)), receiveBuffers(PPM->Size(commSplit)) {}

ExchangeBuffer::~ExchangeBuffer() { Destruct(); }

void ExchangeBuffer::Destruct() {
  ClearBuffers();
  sendBuffers.clear();
  receiveBuffers.clear();
  numOfMessages.clear();
  destinations.clear();
  sendSizes.clear();
}

bool ExchangeBuffer::IsInitialized() const { return init; }

void ExchangeBuffer::CommunicateSize() {
  for (int q = 0; q < PPM->Size(commSplit); ++q)
    numOfMessages[q] = 0;

  for (int q = 0; q < PPM->Size(commSplit); ++q)
    if (sendSizes[q]) ++(numOfMessages[PPM->Proc(commSplit)]);
  PPM->SumOnCommSplit(numOfMessages, commSplit);

  for (int q = PPM->Size(commSplit); q > 0; --q)
    numOfMessages[q] = numOfMessages[q - 1];

  numOfMessages[0] = 0;

  for (int q = 0; q < PPM->Size(commSplit); ++q)
    numOfMessages[q + 1] += numOfMessages[q];

  destinations.resize(numOfMessages[PPM->Size(commSplit)], 0);
  std::vector<size_t> tmp(numOfMessages[PPM->Size(commSplit)], 0);
  int k = numOfMessages[PPM->Proc(commSplit)];

  for (int q = 0; q < PPM->Size(commSplit); ++q) {
    if (sendSizes[q]) {
      tmp[k] = sendSizes[q];
      destinations[k++] = q;
    }
  }

  sendSizes = tmp;
  PPM->SumOnCommSplit(sendSizes, commSplit);
  PPM->SumOnCommSplit(destinations, commSplit);
  numberOfSendMessages =
      numOfMessages[PPM->Proc(commSplit) + 1] - numOfMessages[PPM->Proc(commSplit)];
  numberOfReceiveMessages = 0;

  for (int i = 0; i < numOfMessages[PPM->Size(commSplit)]; ++i)
    if (destinations[i] == PPM->Proc(commSplit)) numberOfReceiveMessages++;

  init = true;
}

size_t ExchangeBuffer::ReceiveSize(int q) const {
  for (int k = numOfMessages[q]; k < numOfMessages[q + 1]; ++k)
    if (destinations[k] == PPM->Proc(commSplit)) return sendSizes[k];
  return 0;
}

void ExchangeBuffer::SendSize(size_t m, int q) { sendSizes[q] = m; }

int ExchangeBuffer::SendMessages() const { return numberOfSendMessages; }

int ExchangeBuffer::RecvMessages() const { return numberOfReceiveMessages; }

int ExchangeBuffer::Messages(int q) const { return numOfMessages[q]; }

int ExchangeBuffer::MessageDest(int q) const { return destinations[q]; }

size_t ExchangeBuffer::MessageSize(int q) const { return sendSizes[q]; }

void ExchangeBuffer::CommunicateSizes() {
  for (short q = 0; q < PPM->Size(commSplit); ++q)
    SendSize(Send(q).size(), q);
  CommunicateSize();
  for (short q = 0; q < PPM->Size(commSplit); ++q)
    Receive(q).resize(ReceiveSize(q));
}

Buffer &ExchangeBuffer::Send(int q) { return sendBuffers[q]; }

Buffer &ExchangeBuffer::Receive(int q) { return receiveBuffers[q]; }

void ExchangeBuffer::Communicate() {
  if (!IsInitialized()) CommunicateSizes();
  PPM->Communicate(*this, commSplit);
  for (int i = 0; i < sendBuffers.size(); ++i)
    sendBuffers[i].Destruct();
}

ExchangeBuffer &ExchangeBuffer::Rewind() {
  for (short q = 0; q < PPM->Size(commSplit); ++q) {
    Send(q).rewind();
    Receive(q).rewind();
  }
  return *this;
}

void ExchangeBuffer::ClearBuffers() {
  for (int i = 0; i < sendBuffers.size(); ++i)
    sendBuffers[i].Destruct();
  for (int i = 0; i < receiveBuffers.size(); ++i)
    receiveBuffers[i].Destruct();
}

int ExchangeBuffer::CommSplit() const { return commSplit; }

std::ostream &operator<<(std::ostream &os, const ExchangeBuffer &exBuffer) {
  os << " on " << PPM->Proc(exBuffer.CommSplit()) << " :"
     << "send " << exBuffer.SendMessages() << " : ";

  for (int k = exBuffer.Messages(PPM->Proc(exBuffer.CommSplit()));
       k < exBuffer.Messages(PPM->Proc(exBuffer.CommSplit()) + 1); ++k)
    os << exBuffer.MessageDest(k) << "|" << exBuffer.MessageSize(k) << " ";

  os << "recv " << exBuffer.RecvMessages() << " : ";

  for (int q = 0; q < PPM->Size(exBuffer.CommSplit()); ++q)
    for (int k = exBuffer.Messages(q); k < exBuffer.Messages(q + 1); ++k)
      if (exBuffer.MessageDest(k) == PPM->Proc(exBuffer.CommSplit()))
        os << q << "|" << exBuffer.MessageSize(k) << " ";

  return os << endl;
}
