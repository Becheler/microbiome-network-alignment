#ifndef ORCA_H
#define ORCA_H

#include <string>
#include <algorithm> // std::min, std::max
#include <unordered_map> // std::unordered_map

///
/// @brief Encapsulate Orca computation logic
///
class ORCA
{
  ///
  /// @brief Nested class for a hashable key for ordered pairs
  ///
  class key_pair
  {
  private:
    // type alias
    using self_type = key_pair;
    // data memebers
    int a;
    int b;
  public:
    ///
    /// @brief Constructor using initialization list
    ///
    key_pair(int a0, int b0): a( std::min(a0,b0) ), b( std::max(a0,b0) ){}
    ///
    /// @brief Comparison operator
    ///
    bool operator<(const self_type &other) const
    {
      if (this->a == other.a)
      return this->b < other.b;
      else
      return this->a < other.a;
    }
    ///
    /// @brief Equality Comparable
    ///
    bool operator==(const self_type &other) const
    {
      return this->a == other.a && this->b == other.b;
    }
    ///
    /// @brief Hashable and usable in hashmaps
    ///
    struct hash
    {
      size_t operator()(const key_pair &x) const
      {
        return (x.a << 8) ^ (x.b << 0);
      }
    };
  }; // end class key_pair

  ///
  /// @brief Nested class for a hashable key for ordered triple
  ///
  class key_triple
  {
  private:
    // type alias
    using self_type = key_triple;
    // data members
    int a;
    int b;
    int c;
  public:
    ///
    /// @brief Constructor using initialization list
    ///
    key_triple(int a0, int b0, int c0):
    a(a0),
    b(b0),
    c(c0)
    {
      if (a > b) swap(a, b);
      if (b > c) swap(b, c);
      if (a > b) swap(a, b);
    }
    ///
    /// @brief Comparison operator
    ///
    bool operator<(const self_type &other) const
    {
      if (this->a == other.a)
      {
        if (this->b == other.b)
        return this->c < other.c;
        else
        return this->b < other.b;
      } else {
        return this->a < other.a;
      }
    }
    ///
    /// @brief EqualityComparable
    ///
    bool operator==(const self_type &other) const
    {
      return this->a == other.a &&
      this->b == other.b &&
      this->c == other.c;
    }
    ///
    /// @brief Hashable and usable in hashmaps
    ///
    struct hash
    {
      size_t operator()(const key_triple &x) const
      {
        return (x.a << 16) ^ (x.b << 8) ^ (x.c << 0);
      }
    };
  }; // end class key_triple

  // Type aliases for readibility
  using common2_type = std::unordered_map<key_pair, int, key_pair::hash>;
  using common3_type = std::unordered_map<key_triple, int, key_triple::hash>;

  // ORCA alorithm state (member data)
  common2_type common2;
  common3_type common3;
  
  common2_type::iterator common2_it;
  common3_type::iterator common3_it;

}; // end class ORCA
#endif
