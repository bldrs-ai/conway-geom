/* MPL License: https://github.com/nickcastel50/conway-geom/blob/typescript_api/LICENSE.md */

#include "IfcTokenStream.h"

namespace webifc::parsing {

IfcTokenStream::IfcFileStream::IfcFileStream(
    const std::function<uint32_t(char *, size_t, size_t)> &requestData,
    uint32_t size)
    : _dataSource(requestData), _size(size) {
  _buffer = nullptr;
  load();
}

void IfcTokenStream::IfcFileStream::load() {
  if (_buffer == nullptr) _buffer = new char[_size];
  prev = _buffer[_currentSize - 1];
  _currentSize = _dataSource(_buffer, _startRef, _size);
  _pointer = 0;
}

void IfcTokenStream::IfcFileStream::Go(uint32_t ref) {
  _startRef = ref;
  load();
}

void IfcTokenStream::IfcFileStream::Forward() {
  _pointer++;
  if (_pointer == _currentSize && _currentSize != 0) {
    _startRef += _currentSize;
    load();
  }
}

void IfcTokenStream::IfcFileStream::Back() {
  _pointer--;
  if (_pointer < 0) {
    _startRef--;
    load();
  }
}

void IfcTokenStream::IfcFileStream::Clear() {
  delete _buffer;
  _buffer = nullptr;
}

char IfcTokenStream::IfcFileStream::Prev() {
  if (_pointer == 0) return prev;
  return _buffer[_pointer - 1];
}

bool IfcTokenStream::IfcFileStream::IsAtEnd() {
  return _pointer == _currentSize && _currentSize == 0;
}

size_t IfcTokenStream::IfcFileStream::GetRef() { return _startRef + _pointer; }

char IfcTokenStream::IfcFileStream::Get() { return _buffer[_pointer]; }
}  // namespace webifc::parsing
