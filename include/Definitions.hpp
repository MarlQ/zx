#ifndef ZX_INCLUDE_DEFINITIONS_HPP_
#define ZX_INCLUDE_DEFINITIONS_HPP_

#include <stdexcept>

namespace zx {
enum class EdgeType { Simple, Hadamard };
enum class VertexType { Boundary, Z, X };
using Vertex = int32_t;
using Col = int32_t;
using Qubit = std::int_fast8_t;
using fp = double;


constexpr double MAX_DENOM = 1e9; // TODO: maybe too high
constexpr double PARAMETER_TOLERANCE = 1e-13;
constexpr double TOLERANCE = 1e-13;
static constexpr double PI =
    3.141592653589793238462643383279502884197169399375105820974L;

class ZXException : public std::invalid_argument {
  std::string msg;

public:
  explicit ZXException(std::string msg)
      : std::invalid_argument("ZX Exception"), msg(std::move(msg)) {}

  [[nodiscard]] const char *what() const noexcept override {
    return msg.c_str();
  }
};
} // namespace zx
#endif /* JKQZX_INCLUDE_DEFINITIONS_HPP_ */
