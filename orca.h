#ifndef ORCA_H
#define ORCA_H

#include <string>
#include <algorithm> // std::min, std::max, std::binary_search
#include <utility> // std::min, std::max, std::binary_search
#include <unordered_map> // std::unordered_map
#include <fstream> // ofstream
#include <iostream> // std::cout
#include <exception> // std::invalid_argument
#include <functional> // std::function
#include <iterator> // std::advance
#include <cassert>
#include <set>
///
/// @brief Nested class for a hashable key for pairs of nodes
///
class key_pair
{
private:
  // type alias
  using self_type = key_pair;
public:
  // data memebers
  int a;
  int b;
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
/// @brief Nested class for a hashable key for triplet of nodes
///
class key_triple
{
private:
  // type alias
  using self_type = key_triple;
public:
  // data members
  int a;
  int b;
  int c;
  ///
  /// @brief Constructor using initialization list
  ///
  key_triple(int a0, int b0, int c0): a(a0), b(b0), c(c0)
  {
    if (a > b) std::swap(a, b);
    if (b > c) std::swap(b, c);
    if (a > b) std::swap(a, b);
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
    return this->a == other.a && this->b == other.b && this->c == other.c;
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

///
/// @brief Exchangeable implementations of an adjacent matrix (policy-based design)
///
namespace adjacent_policy
{
  // or `long long` or  `boost::multiprecision::cpp_int` (versatile) or `boost::multiprecision::cpp_int` (fast)
  using big_int = unsigned int;
  ///
  /// @brief adj is a vector of vector
  ///
  class default_impl
  {
  public:
    // type used for composition in the ORCA host class
    using adjacent_matrix_type = std::vector<std::vector<int>>;
    // Precision for very large number of nodes
    using big_int = typename adjacent_policy::big_int;
  private:
    // Reference on edges
    const std::vector<int>& _deg;
    // adj[x] - adjacency list of node x
    adjacent_matrix_type _adj;
  public:
    ///
    /// @brief Read only access
    ///
    const adjacent_matrix_type & get_adjacent_matrix() const { return _adj;}
    ///
    /// @brief Policy constructor
    ///
    default_impl(int n, int, const std::vector<key_pair>, const std::vector<int> & deg):
    _deg(deg),
    _adj(build_implementation(n, deg)){}
    ///
    /// @brief Checks if an element equivalent to value appears within the range x, y
    ///
    /// @note common interface between policies
    ///
    bool are_adjacent(int x, int y) const
    {
      // Preconditions
      assert(x < _adj.size());
      assert( _deg.at(x) < _adj.at(x).size() );

      auto first = _adj.at(x).cbegin();
      auto last = first;
      std::advance(last, _deg.at(x));
      return std::binary_search(first, last, y);
    }
    ///
    /// @brief Harmonize information access through different policies
    ///
    auto operator()(int x, int y) const
    {
      // Preconditions
      assert(x < _adj.size() );
      assert(y < _adj.at(x).size() );
      return _adj.at(x).at(y);
    }
    ///
    /// @brief Harmonize information access through different policies
    ///
    auto& operator()(int x, int y)
    {
      // Preconditions
      assert(x < _adj.size());
      assert(y < _adj.at(x).size());
      return _adj.at(x).at(y);
    }
    ///
    /// @brief Sort elements in a range along the 2nd dimension
    ///
    void sort(int i, int j)
    {
      // Preconditions
      assert(i < _adj.size());
      assert(j < _adj.size());
      assert(i + j < _adj.size());
      // iterator on first elements of ith vector second dimension
      auto first = _adj.at(i).begin();
      auto last  = _adj.at(i).begin() + j;
      std::sort(first, last);
    }
  private:
    ///
    /// @brief Build implementation
    ///
    std::vector<std::vector<int>> build_implementation(int n, const std::vector<int> & deg)
    {
      std::vector<std::vector<int>> adj(n);
      for (std::size_t i = 0; auto& it : adj)
      {
        it.resize(deg.at(i));
      }
      return adj;
    }
  };

  ///
  /// @brief adj is one long vector
  ///
  class compressed
  {
  public:
    // type used for composition in the ORCA host class
    using adjacent_matrix_type = std::vector<int>;
    // Precision for very large number of nodes
    using big_int = adjacent_policy::big_int;
  private:
    // number of nodes
    int _n;
    // chunk size
    static constexpr int chunk = 8 * sizeof(int);
    // compressed adjacency matrix
    adjacent_matrix_type _adj;
  public:
    //
    ///
    /// @ brief Build the adjacent matrix
    ///
    /// @note Strategy is set up adjacency matrix if it's smaller than 100MB
    ///
    std::vector<int> build_implementation(int n, int m, const std::vector<key_pair> &edges) const
    {
      // Initialize a zero vector
      std::vector<int> adj(n*n, 0);
      // Set matrix to 1 since a and b are connected
      for (const auto & edge : edges)
      {
        adj.at(edge.a * n + edge.b) = 1;
        adj.at(edge.b * n + edge.a) = 1;
      }
      return adj;
    }
  public:
    ///
    /// @brief Checks if an element equivalent to value appears within the range x, y
    ///
    /// @note common interface between policies
    ///
    bool are_adjacent(int x, int y) const
    {
      return _adj[(x * _n + y) / chunk] & (1 << ((x * _n + y) % chunk));
    }
    ///
    /// @brief Harmonize information access through different policies
    ///
    auto operator()(int x, int y) const
    {
      size_t index = x * _n + y;
      assert(index < _adj.size());
      return _adj.at(index);
    }
    ///
    /// @brief Harmonize information access through different policies
    ///
    auto& operator()(int x, int y)
    {
      size_t index = x * _n + y;
      assert(index < _adj.size());
      return  _adj.at(x * _n + y );
    }
    ///
    /// @brief Sort elements in a range
    ///
    void sort(int i, int j)
    {
      auto index = i * _n ;
      std::cout << i << " " << j << " " << index << " " << _adj.size() << std::endl;
      assert( index < _adj.size());
      auto first = _adj.begin() + index ;
      auto last = _adj.begin() + index + j;
      std::sort(first, last);
    }
    ///
    /// @brief Policy constructor
    ///
    compressed(int n, int m, const std::vector<key_pair> &edges, const std::vector<int> &):
    _n(n),
    _adj(build_implementation(n, m, edges)){}
    ///
    /// @brief Should the user pick this implementation?
    static inline constexpr bool should_use_based_on(int nb_nodes)
    {
      return static_cast<big_int>(nb_nodes * nb_nodes) < 100LL * 1024 * 1024 * 8;
    }

  };
} // end namespace adjacent_implementation

///
/// @brief Orbit Counting Algorithm
///
/// @note Encapsulate Orca computation logic
///
template<class Strategy=adjacent_policy::default_impl>
class ORCA
{
private:
  // Type aliases for readibility
  using common2_type = std::unordered_map<key_pair, int, key_pair::hash>;
  using common3_type = std::unordered_map<key_triple, int, key_triple::hash>;
  using strategy_type = Strategy;
  using big_int = strategy_type::big_int;
  // number of nodes
  int _n;
  // number of edges
  int _m;
  // degrees of individual nodes
  std::vector<int> _deg;
  // list of edge
  std::vector<key_pair> _edges;
  // Implementation of adjacent matrix is variable
  strategy_type _policy;
  // inc[x] - incidence list of node x: (y, edge id)
  std::vector<std::vector<std::pair<int,int>>> _inc;
  // orbit[x][o] - how many times does node x participate in orbit o
  std::vector<std::vector<big_int>> _orbits;
  // ORCA alorithm state (member data)
  common2_type _common2;
  // stores the number of nodes that are adjacent to some triplets of nodes
  common3_type _common3;
  // stores the triangles graphlets
  std::vector<int> _triangles;
  // stores graphlets
  std::vector<big_int> _C5;

  ///
  /// @brief Return the value stored at the key, else 0
  ///
  /// @note avoid inflation of the map due to lookups of absent elements.
  ///
  template<class key, class T>
  T find_or_zero(key x, const std::unordered_map<key, T, typename key::hash>& map) const
  {
    auto it = map.find(x);
    return (it == map.end()) ? it->second : 0;
  }

  ///
  /// @brief Number of nodes that are adjacent to some pair of nodes
  ///
  int common2_get(const key_pair & x) const
  {
    return find_or_zero(x, this->_common2);
  }

  ///
  /// @brief Number of nodes that are adjacent to some triplets of nodes
  ///
  int common3_get(const key_triple & x) const
  {
    return find_or_zero(x, this->_common3);
  }

  ///
  /// @brief Initialize an incidence matrix
  ///
  auto resize_incidence_matrix(int n) const
  {
    std::vector< std::vector< std::pair<int,int> >> inc(n);
    for(auto &it : inc)
    {
      it.resize(n);
    }
    return inc;
  }

  ///
  /// @brief Initialize a the orbit atrix
  ///
  auto resize_orbits(big_int n) const
  {
    std::vector<std::vector<big_int>> orbits(n);
    for (auto & it : orbits)
    {
      it.resize(73);
    }
    return orbits;
  }

  ///
  /// @brief Sort elements in a specific range [first, last]
  ///
  void sort_incidence_matrix_range(int i, int j)
  {
    // Preconditions
    auto first = _inc.at(i).begin();
    auto last =  _inc.at(i).begin() + j;
    std::sort(first, last);
  }
public:
  /// @brief Move constructor
  ///
  /// @note Move semantics: edges and deg content are "stolen" from the previous context
  ///
  ORCA(int n, int m, std::vector<key_pair> &&edges, std::vector<int> &&deg) noexcept :
  _n(n),
  _m(m),
  _deg(std::move(deg)),
  _edges(std::move(edges)),
  _policy(n, m, _edges, _deg),
  _inc(resize_incidence_matrix(n)),
  _orbits(resize_orbits(n)),
  _triangles(m),
  _C5(n)
  {
    std::vector<int> d(n);

    for (std::size_t i = 0; auto const & it : edges)
    {
      int a = it.a;
      int b = it.b;

      // accessing policy protected data, not great but okay-ish
      this->_policy(a, d[a]) = b;
      this->_policy(b, d[b]) = a;

      this->_inc[a][d[a]] = std::pair<int,int>(b, i);
      this->_inc[b][d[b]] = std::pair<int,int>(a, i);

      d[a]++;
      d[b]++;
    }

    for (int i = 0; i < _n; i++)
    {
      this->_policy.sort(i, _deg[i]);
      sort_incidence_matrix_range(i, _deg[i]);
    }
  } // end constructor
  ///
  /// @brief Write accumulated results to a file
  ///
  void write_results_to(const std::string &output_file) const
  {
    std::ofstream myfile;
    myfile.open (output_file);
    myfile << format();
    myfile.close();
  }
  ///
  /// @brief Precompute common nodes
  ///
  void precompute_common_nodes()
  {
    // alias on policy operator (lambda function)
    auto are_adjacent = [this](int x, int y){ return this->_policy.are_adjacent(x, y);};

    // precompute common nodes
    std::cout << "stage 1 - precomputing common nodes" << std::endl ;
    int frac_prev = -1;
    for (int x = 0; x < _n; x++)
    {
      auto frac = 100LL * x / _n;

      if (frac != frac_prev)
      {
        std::cout << frac << std::endl;
        frac_prev = frac;
      }

      for (int n1 = 0; n1 < _deg.at(x); n1++)
      {
        int a = _policy(x, n1);
        for (int n2 = n1 + 1; n2 < _deg.at(x); n2++)
        {
          int b = _policy(x, n2);
          auto ab = key_pair(a, b);
          this->_common2[ab]++;
          for (int n3 = n2 + 1; n3 < _deg[x]; n3++)
          {
            int c = _policy(x, n3);
            int st = are_adjacent(a, b) + are_adjacent(a, c) + are_adjacent(b, c);
            if (st < 2) continue;
            auto abc = key_triple(a, b, c);
            this->_common3[abc]++;
          }
        }
      }
    }
  }
  ///
  /// @brief Precompute triangles that span over edges
  ///
  void precompute_triangles_that_span_over_edges()
  {
    for (int i = 0; i < _m; i++)
    {
      int x = _edges.at(i).a;
      int y = _edges.at(i).b;
      for (int xi = 0, yi = 0; xi < _deg.at(x) && yi < _deg.at(y);)
      {
        if (_policy(x, xi) == _policy(y, yi))
        {
          _triangles.at(i)++;
          xi++;
          yi++;
        } else if ( _policy(x,xi) < _policy(y,yi) ) {
          xi++;
        } else {
          yi++;
        }
      }
    }
  }
  ///
  /// @brief Count full graphlets
  ///
  void count_full_graphlets(int n)
  {
    std::cout << "stage 2 - counting full graphlets" << std::endl;

    std::vector<int> neigh(n);
    std::vector<int> neigh2(n);

    // alias on policy operator (lambda function)
    auto are_adjacent = [this](int x, int y){ return _policy.are_adjacent(x, y);};

    int nn;
    int nn2;
    int frac_prev = -1;

    for (int x = 0; x < n; x++)
    {
      auto frac = 100LL * x / n;
      if (frac != frac_prev)
      {
        std::cout << frac << std::endl;
        frac_prev = frac;
      }
      for (int nx = 0; nx < _deg.at(x); nx++)
      {
        int y = _policy(x, nx);
        if (y >= x) break;
        nn = 0;
        for (int ny = 0; ny < _deg.at(y); ny++)
        {
          int z = _policy(y, ny);
          if (z >= y) break;
          if (are_adjacent(x, z))
          {
            neigh.at(nn++) = z;
          }
        }
        for (int i = 0; i < nn; i++)
        {
          int z = neigh.at(i);
          nn2 = 0;
          for (int j = i + 1; j < nn; j++)
          {
            int zz = neigh.at(j);
            if (are_adjacent(z, zz))
            {
              neigh2.at(nn2++) = zz;
            }
          }
          for (int i2 = 0; i2 < nn2; i2++)
          {
            int zz = neigh2.at(i2);
            for (int j2 = i2 + 1; j2 < nn2; j2++)
            {
              int zzz = neigh2.at(j2);
              if (are_adjacent(zz, zzz))
              {
                _C5.at(x)++;
                _C5.at(y)++;
                _C5.at(z)++;
                _C5.at(zz)++;
                _C5.at(zzz)++;
              }
            }
          }
        }
      }
    }
  }
  ///
  /// @brief Set up a system of equations relating orbit counts
  ///
  // void set_up_equations()
  // {
  //
  // }
  /// @brief Count graphlets on max 5 nodes
  ///
  void count_orbits()
  {

    clock_t startTime, endTime;
    startTime = clock();

    clock_t startTime_all, endTime_all;
    startTime_all = startTime;

    precompute_common_nodes();
    precompute_triangles_that_span_over_edges();

    endTime = clock();
    std::cout << (double)(endTime - startTime) / CLOCKS_PER_SEC << " sec" << std::endl;
    startTime = endTime;

    count_full_graphlets(_n);

    // set up a system of equations relating orbit counts
    // solve equations

    endTime = clock();
    std::cout << (double)(endTime - startTime) / CLOCKS_PER_SEC << " sec" << std::endl;
    startTime = endTime;

    std::vector<int> common_x(_n);
    std::vector<int> common_x_list(_n);
    int ncx = 0;
    std::vector<int> common_a(_n);
    std::vector<int> common_a_list(_n);
    int nca = 0;

    // set up a system of equations relating orbit counts
    printf("stage 3 - building systems of equations\n");
    int frac_prev = -1;

    for (int x = 0; x < _n; x++)
    {

      auto frac = 100LL * x / _n;
      if (frac != frac_prev)
      {
        std::cout << frac << std::endl;
        frac_prev = frac;
      }

      for (int i = 0; i < ncx; i++) common_x.at(common_x_list.at(i)) = 0;
      ncx = 0;

      // smaller graphlets
      _orbits.at(x).at(0) = _deg.at(x);
      for (int nx1 = 0; nx1 < _deg.at(x); nx1++)
      {
        int a = _policy(x, nx1);

        for (int nx2 = nx1 + 1; nx2 < _deg.at(x); nx2++)
        {
          int b = _policy(x, nx2);
          if (_policy.are_adjacent(a, b))
            _orbits.at(x).at(3)++;
          else
            _orbits.at(x).at(2)++;
        }

        for (int na = 0; na < _deg.at(a); na++)
        {
          int b = _policy(a, na);
          if (b != x && !_policy.are_adjacent(x, b))
          {
            _orbits.at(x).at(1)++;
            if (common_x.at(b) == 0) common_x_list.at(ncx++) = b;
            common_x.at(b)++;
          }
        }
      }

      big_int f_71 = 0, f_70 = 0, f_67 = 0, f_66 = 0, f_58 = 0, f_57 = 0;                                // 14
      big_int f_69 = 0, f_68 = 0, f_64 = 0, f_61 = 0, f_60 = 0, f_55 = 0, f_48 = 0, f_42 = 0, f_41 = 0;  // 13
      big_int f_65 = 0, f_63 = 0, f_59 = 0, f_54 = 0, f_47 = 0, f_46 = 0, f_40 = 0;                      // 12
      big_int f_62 = 0, f_53 = 0, f_51 = 0, f_50 = 0, f_49 = 0, f_38 = 0, f_37 = 0, f_36 = 0;            // 8
      big_int f_44 = 0, f_33 = 0, f_30 = 0, f_26 = 0;                                                    // 11
      big_int f_52 = 0, f_43 = 0, f_32 = 0, f_29 = 0, f_25 = 0;                                          // 10
      big_int f_56 = 0, f_45 = 0, f_39 = 0, f_31 = 0, f_28 = 0, f_24 = 0;                                // 9
      big_int f_35 = 0, f_34 = 0, f_27 = 0, f_18 = 0, f_16 = 0, f_15 = 0;                                // 4
      big_int f_17 = 0;                                                                                  // 5
      big_int f_22 = 0, f_20 = 0, f_19 = 0;                                                              // 6
      big_int f_23 = 0, f_21 = 0;                                                                        // 7

      for (int nx1 = 0; nx1 < _deg.at(x); nx1++)
      {
        auto a = _inc.at(x).at(nx1).first;
        auto xa = _inc.at(x).at(nx1).second;

        for (int i = 0; i < nca; i++) common_a.at(common_a_list.at(i)) = 0;
        nca = 0;

        for (int na = 0; na < _deg.at(a); na++)
        {
          auto b = _policy(a, na);
          for (int nb = 0; nb < _deg.at(b); nb++)
          {
            int c = _policy(b, nb);

            if (c == a || _policy.are_adjacent(a, c)) continue;

            if (common_a.at(c) == 0) common_a_list.at(nca++) = c;
            common_a.at(c)++;
          }
        }

        // x = orbit-14 (tetrahedron)
        for (int nx2 = nx1 + 1; nx2 < _deg.at(x); nx2++)
        {
          auto b = _inc.at(x).at(nx2).first;
          auto xb = _inc[x][nx2].second;

          if (!_policy.are_adjacent(a, b)) continue;

          for (int nx3 = nx2 + 1; nx3 < _deg.at(x); nx3++)
          {
            auto c = _inc.at(x).at(nx3).first;
            auto xc = _inc.at(x).at(nx3).second;

            if (!_policy.are_adjacent(a, c) || !_policy.are_adjacent(b, c)) continue;

            _orbits.at(x).at(14)++;
            f_70 += common3_get(key_triple(a, b, c)) - 1;

            // // debug
            // if (nx2 == nx1 + 6)
            //     printf("%d %d %d\n", x, a, b);

            f_71 += (_triangles.at(xa) > 2 && _triangles.at(xb) > 2) ? (common3_get(key_triple(x, a, b)) - 1) : 0;
            f_71 += (_triangles.at(xa) > 2 && _triangles.at(xc) > 2) ? (common3_get(key_triple(x, a, c)) - 1) : 0;
            f_71 += (_triangles.at(xb) > 2 && _triangles.at(xc) > 2) ? (common3_get(key_triple(x, b, c)) - 1) : 0;
            f_67 += _triangles.at(xa) - 2 + _triangles.at(xb) - 2 + _triangles.at(xc) - 2;
            f_66 += common2_get(key_pair(a, b)) - 2;
            f_66 += common2_get(key_pair(a, c)) - 2;
            f_66 += common2_get(key_pair(b, c)) - 2;
            f_58 += _deg.at(x) - 3;
            f_57 += _deg.at(a) - 3 + _deg.at(b) - 3 + _deg.at(c) - 3;
          }
        }

        // x = orbit-13 (diamond)
        for (int nx2 = 0; nx2 < _deg.at(x); nx2++)
        {
          auto b = _inc.at(x).at(nx2).first;
          auto xb = _inc.at(x).at(nx2).second;

          if (!_policy.are_adjacent(a, b)) continue;

          for (int nx3 = nx2 + 1; nx3 < _deg.at(x); nx3++)
          {
            auto c = _inc.at(x).at(nx3).first;
            auto xc = _inc.at(x).at(nx3).second;

            if (!_policy.are_adjacent(a, c) || _policy.are_adjacent(b, c)) continue;

            _orbits.at(x).at(13)++;
            f_69 += (_triangles.at(xb) > 1 && _triangles.at(xc) > 1) ? (common3_get(key_triple(x, b, c)) - 1) : 0;
            f_68 += common3_get(key_triple(a, b, c)) - 1;
            f_64 += common2_get(key_pair(b, c)) - 2;
            f_61 += _triangles.at(xb) - 1 + _triangles.at(xc) - 1;
            f_60 += common2_get(key_pair(a, b)) - 1;
            f_60 += common2_get(key_pair(a, c)) - 1;
            f_55 += _triangles.at(xa) - 2;
            f_48 += _deg.at(b) - 2 + _deg.at(c) - 2;
            f_42 += _deg.at(x) - 3;
            f_41 += _deg.at(a) - 3;
          }
        }

        // x = orbit-12 (diamond)
        for (int nx2 = nx1 + 1; nx2 < _deg.at(x); nx2++)
        {
          int b = _inc.at(x).at(nx2).first;

          if (!_policy.are_adjacent(a, b)) continue;

          for (int na = 0; na < _deg.at(a); na++)
          {
            auto c = _inc.at(a).at(na).first;
            auto ac =_inc.at(a).at(na).second;

            if (c == x || _policy.are_adjacent(x, c) || !_policy.are_adjacent(b, c)) continue;

            _orbits.at(x).at(12)++;
            f_65 += (_triangles.at(ac) > 1) ? common3_get(key_triple(a, b, c)) : 0;
            f_63 += common_x.at(c) - 2;
            f_59 += _triangles.at(ac) - 1 + common2_get(key_pair(b, c)) - 1;
            f_54 += common2_get(key_pair(a, b)) - 2;
            f_47 += _deg.at(x) - 2;
            f_46 += _deg.at(c) - 2;
            f_40 += _deg.at(a) - 3 + _deg.at(b) - 3;
          }
        }

        // x = orbit-8 (cycle)
        for (int nx2 = nx1 + 1; nx2 < _deg.at(x); nx2++)
        {
          auto b = _inc.at(x).at(nx2).first;
          auto xb = _inc.at(x).at(nx2).second;

          if (_policy.are_adjacent(a, b)) continue;

          for (int na = 0; na < _deg.at(a); na++)
          {
            auto c = _inc.at(a).at(na).first;
            auto ac = _inc.at(a).at(na).second;

            if (c == x || _policy.are_adjacent(x, c) || !_policy.are_adjacent(b, c)) continue;

            _orbits.at(x).at(8)++;
            f_62 += (_triangles.at(ac) > 0) ? common3_get(key_triple(a, b, c)) : 0;
            f_53 += _triangles.at(xa) + _triangles.at(xb);
            f_51 += _triangles.at(ac) + common2_get(key_pair(c, b));
            f_50 += common_x.at(c) - 2;
            f_49 += common_a.at(b) - 2;
            f_38 += _deg.at(x) - 2;
            f_37 += _deg.at(a) - 2 + _deg.at(b) - 2;
            f_36 += _deg.at(c) - 2;
          }
        }

        // x = orbit-11 (paw)
        for (int nx2 = nx1 + 1; nx2 < _deg.at(x); nx2++)
        {
          int b = _inc.at(x).at(nx2).first;

          if (!_policy.are_adjacent(a, b)) continue;

          for (int nx3 = 0; nx3 < _deg.at(x); nx3++)
          {
            auto c = _inc.at(x).at(nx3).first;
            auto xc = _inc.at(x).at(nx3).second;

            if (c == a || c == b || _policy.are_adjacent(a, c) || _policy.are_adjacent(b, c)) continue;

            _orbits.at(x).at(11)++;
            f_44 += _triangles.at(xc);
            f_33 += _deg.at(x) - 3;
            f_30 += _deg.at(c) - 1;
            f_26 += _deg.at(a) - 2 + _deg.at(b) - 2;
          }
        }

        // x = orbit-10 (paw)
        for (int nx2 = 0; nx2 < _deg.at(x); nx2++)
        {
          int b = _inc.at(x).at(nx2).first;

          if (!_policy.are_adjacent(a, b)) continue;

          for (int nb = 0; nb < _deg.at(b); nb++)
          {
            auto c = _inc.at(b).at(nb).first;
            auto bc = _inc.at(b).at(nb).second;

            if (c == x || c == a || _policy.are_adjacent(a, c) || _policy.are_adjacent(x, c)) continue;

            _orbits.at(x).at(10)++;
            f_52 += common_a.at(c) - 1;
            f_43 += _triangles.at(bc);
            f_32 += _deg.at(b) - 3;
            f_29 += _deg.at(c) - 1;
            f_25 += _deg.at(a) - 2;
          }
        }

        // x = orbit-9 (paw)
        for (int na1 = 0; na1 < _deg.at(a); na1++)
        {
          auto b = _inc.at(a).at(na1).first;
          auto ab = _inc.at(a).at(na1).second;

          if (b == x || _policy.are_adjacent(x, b)) continue;

          for (int na2 = na1 + 1; na2 < _deg.at(a); na2++)
          {
            auto c = _inc.at(a).at(na2).first;
            auto ac = _inc.at(a).at(na2).second;

            if (c == x || !_policy.are_adjacent(b, c) || _policy.are_adjacent(x, c)) continue;

            _orbits.at(x).at(9)++;
            f_56 += (_triangles.at(ab) > 1 && _triangles.at(ac) > 1) ? common3_get(key_triple(a, b, c)) : 0;
            f_45 += common2_get(key_pair(b, c)) - 1;
            f_39 += _triangles.at(ab) - 1 + _triangles.at(ac) - 1;
            f_31 += _deg.at(a) - 3;
            f_28 += _deg.at(x) - 1;
            f_24 += _deg.at(b) - 2 + _deg.at(c) - 2;
          }
        }

        // x = orbit-4 (path)
        for (int na = 0; na < _deg.at(a); na++)
        {
          auto b = _inc.at(a).at(na).first;

          if (b == x || _policy.are_adjacent(x, b)) continue;

          for (int nb = 0; nb < _deg.at(b); nb++)
          {
            int c = _inc.at(b).at(nb).first;
            auto bc = _inc.at(b).at(nb).second;

            if (c == a || _policy.are_adjacent(a, c) || _policy.are_adjacent(x, c)) continue;

            _orbits.at(x).at(4)++;
            f_35 += common_a.at(c) - 1;
            f_34 += common_x.at(c);
            f_27 += _triangles.at(bc);
            f_18 += _deg.at(b) - 2;
            f_16 += _deg.at(x) - 1;
            f_15 += _deg.at(c) - 1;
          }
        }

        // x = orbit-5 (path)
        for (int nx2 = 0; nx2 < _deg.at(x); nx2++)
        {
          int b = _inc.at(x).at(nx2).first;

          if (b == a || _policy.are_adjacent(a, b)) continue;

          for (int nb = 0; nb < _deg.at(b); nb++)
          {
            int c = _inc.at(b).at(nb).first;

            if (c == x || _policy.are_adjacent(a, c) || _policy.are_adjacent(x, c)) continue;

            _orbits.at(x).at(5)++;
            f_17 += _deg.at(a) - 1;
          }
        }

        // x = orbit-6 (claw)
        for (int na1 = 0; na1 < _deg.at(a); na1++)
        {
          int b = _inc.at(a).at(na1).first;

          if (b == x || _policy.are_adjacent(x, b)) continue;

          for (int na2 = na1 + 1; na2 < _deg.at(a); na2++)
          {
            int c = _inc.at(a).at(na2).first;
            if (c == x || _policy.are_adjacent(x, c) || _policy.are_adjacent(b, c)) continue;
            _orbits.at(x).at(6)++;
            f_22 += _deg.at(a) - 3;
            f_20 += _deg.at(x) - 1;
            f_19 += _deg.at(b) - 1 + _deg.at(c) - 1;
          }
        }

        // x = orbit-7 (claw)
        for (int nx2 = nx1 + 1; nx2 < _deg.at(x); nx2++)
        {
          int b = _inc.at(x).at(nx2).first;
          if (_policy.are_adjacent(a, b)) continue;
          for (int nx3 = nx2 + 1; nx3 < _deg.at(x); nx3++)
          {
            int c = _inc.at(x).at(nx3).first;
            if (_policy.are_adjacent(a, c) || _policy.are_adjacent(b, c)) continue;
            _orbits.at(x).at(7)++;
            f_23 += _deg.at(x) - 3;
            f_21 += _deg.at(a) - 1 + _deg.at(b) - 1 + _deg.at(c) - 1;
          }
        }
      }

      // read and write data accumulated by the algorithm
      auto& orbit = this->_orbits;
      // solve equations
      orbit[x][72] = _C5[x];
      orbit[x][71] = (f_71 - 12 * orbit[x][72]) / 2;
      orbit[x][70] = (f_70 - 4 * orbit[x][72]);
      orbit[x][69] = (f_69 - 2 * orbit[x][71]) / 4;
      orbit[x][68] = (f_68 - 2 * orbit[x][71]);
      orbit[x][67] = (f_67 - 12 * orbit[x][72] - 4 * orbit[x][71]);
      orbit[x][66] = (f_66 - 12 * orbit[x][72] - 2 * orbit[x][71] - 3 * orbit[x][70]);
      orbit[x][65] = (f_65 - 3 * orbit[x][70]) / 2;
      orbit[x][64] = (f_64 - 2 * orbit[x][71] - 4 * orbit[x][69] - 1 * orbit[x][68]);
      orbit[x][63] = (f_63 - 3 * orbit[x][70] - 2 * orbit[x][68]);
      orbit[x][62] = (f_62 - 1 * orbit[x][68]) / 2;
      orbit[x][61] = (f_61 - 4 * orbit[x][71] - 8 * orbit[x][69] - 2 * orbit[x][67]) / 2;
      orbit[x][60] = (f_60 - 4 * orbit[x][71] - 2 * orbit[x][68] - 2 * orbit[x][67]);
      orbit[x][59] = (f_59 - 6 * orbit[x][70] - 2 * orbit[x][68] - 4 * orbit[x][65]);
      orbit[x][58] = (f_58 - 4 * orbit[x][72] - 2 * orbit[x][71] - 1 * orbit[x][67]);
      orbit[x][57] = (f_57 - 12 * orbit[x][72] - 4 * orbit[x][71] - 3 * orbit[x][70] - 1 * orbit[x][67] - 2 * orbit[x][66]);
      orbit[x][56] = (f_56 - 2 * orbit[x][65]) / 3;
      orbit[x][55] = (f_55 - 2 * orbit[x][71] - 2 * orbit[x][67]) / 3;
      orbit[x][54] = (f_54 - 3 * orbit[x][70] - 1 * orbit[x][66] - 2 * orbit[x][65]) / 2;
      orbit[x][53] = (f_53 - 2 * orbit[x][68] - 2 * orbit[x][64] - 2 * orbit[x][63]);
      orbit[x][52] = (f_52 - 2 * orbit[x][66] - 2 * orbit[x][64] - 1 * orbit[x][59]) / 2;
      orbit[x][51] = (f_51 - 2 * orbit[x][68] - 2 * orbit[x][63] - 4 * orbit[x][62]);
      orbit[x][50] = (f_50 - 1 * orbit[x][68] - 2 * orbit[x][63]) / 3;
      orbit[x][49] = (f_49 - 1 * orbit[x][68] - 1 * orbit[x][64] - 2 * orbit[x][62]) / 2;
      orbit[x][48] = (f_48 - 4 * orbit[x][71] - 8 * orbit[x][69] - 2 * orbit[x][68] - 2 * orbit[x][67] - 2 * orbit[x][64] - 2 * orbit[x][61] - 1 * orbit[x][60]);
      orbit[x][47] = (f_47 - 3 * orbit[x][70] - 2 * orbit[x][68] - 1 * orbit[x][66] - 1 * orbit[x][63] - 1 * orbit[x][60]);
      orbit[x][46] = (f_46 - 3 * orbit[x][70] - 2 * orbit[x][68] - 2 * orbit[x][65] - 1 * orbit[x][63] - 1 * orbit[x][59]);
      orbit[x][45] = (f_45 - 2 * orbit[x][65] - 2 * orbit[x][62] - 3 * orbit[x][56]);
      orbit[x][44] = (f_44 - 1 * orbit[x][67] - 2 * orbit[x][61]) / 4;
      orbit[x][43] = (f_43 - 2 * orbit[x][66] - 1 * orbit[x][60] - 1 * orbit[x][59]) / 2;
      orbit[x][42] = (f_42 - 2 * orbit[x][71] - 4 * orbit[x][69] - 2 * orbit[x][67] - 2 * orbit[x][61] - 3 * orbit[x][55]);
      orbit[x][41] = (f_41 - 2 * orbit[x][71] - 1 * orbit[x][68] - 2 * orbit[x][67] - 1 * orbit[x][60] - 3 * orbit[x][55]);
      orbit[x][40] = (f_40 - 6 * orbit[x][70] - 2 * orbit[x][68] - 2 * orbit[x][66] - 4 * orbit[x][65] - 1 * orbit[x][60] - 1 * orbit[x][59] - 4 * orbit[x][54]);
      orbit[x][39] = (f_39 - 4 * orbit[x][65] - 1 * orbit[x][59] - 6 * orbit[x][56]) / 2;
      orbit[x][38] = (f_38 - 1 * orbit[x][68] - 1 * orbit[x][64] - 2 * orbit[x][63] - 1 * orbit[x][53] - 3 * orbit[x][50]);
      orbit[x][37] = (f_37 - 2 * orbit[x][68] - 2 * orbit[x][64] - 2 * orbit[x][63] - 4 * orbit[x][62] - 1 * orbit[x][53] - 1 * orbit[x][51] - 4 * orbit[x][49]);
      orbit[x][36] = (f_36 - 1 * orbit[x][68] - 2 * orbit[x][63] - 2 * orbit[x][62] - 1 * orbit[x][51] - 3 * orbit[x][50]);
      orbit[x][35] = (f_35 - 1 * orbit[x][59] - 2 * orbit[x][52] - 2 * orbit[x][45]) / 2;
      orbit[x][34] = (f_34 - 1 * orbit[x][59] - 2 * orbit[x][52] - 1 * orbit[x][51]) / 2;
      orbit[x][33] = (f_33 - 1 * orbit[x][67] - 2 * orbit[x][61] - 3 * orbit[x][58] - 4 * orbit[x][44] - 2 * orbit[x][42]) / 2;
      orbit[x][32] = (f_32 - 2 * orbit[x][66] - 1 * orbit[x][60] - 1 * orbit[x][59] - 2 * orbit[x][57] - 2 * orbit[x][43] - 2 * orbit[x][41] - 1 * orbit[x][40]) / 2;
      orbit[x][31] = (f_31 - 2 * orbit[x][65] - 1 * orbit[x][59] - 3 * orbit[x][56] - 1 * orbit[x][43] - 2 * orbit[x][39]);
      orbit[x][30] = (f_30 - 1 * orbit[x][67] - 1 * orbit[x][63] - 2 * orbit[x][61] - 1 * orbit[x][53] - 4 * orbit[x][44]);
      orbit[x][29] = (f_29 - 2 * orbit[x][66] - 2 * orbit[x][64] - 1 * orbit[x][60] - 1 * orbit[x][59] - 1 * orbit[x][53] - 2 * orbit[x][52] - 2 * orbit[x][43]);
      orbit[x][28] = (f_28 - 2 * orbit[x][65] - 2 * orbit[x][62] - 1 * orbit[x][59] - 1 * orbit[x][51] - 1 * orbit[x][43]);
      orbit[x][27] = (f_27 - 1 * orbit[x][59] - 1 * orbit[x][51] - 2 * orbit[x][45]) / 2;
      orbit[x][26] = (f_26 - 2 * orbit[x][67] - 2 * orbit[x][63] - 2 * orbit[x][61] - 6 * orbit[x][58] - 1 * orbit[x][53] - 2 * orbit[x][47] - 2 * orbit[x][42]);
      orbit[x][25] = (f_25 - 2 * orbit[x][66] - 2 * orbit[x][64] - 1 * orbit[x][59] - 2 * orbit[x][57] - 2 * orbit[x][52] - 1 * orbit[x][48] - 1 * orbit[x][40]) / 2;
      orbit[x][24] = (f_24 - 4 * orbit[x][65] - 4 * orbit[x][62] - 1 * orbit[x][59] - 6 * orbit[x][56] - 1 * orbit[x][51] - 2 * orbit[x][45] - 2 * orbit[x][39]);
      orbit[x][23] = (f_23 - 1 * orbit[x][55] - 1 * orbit[x][42] - 2 * orbit[x][33]) / 4;
      orbit[x][22] = (f_22 - 2 * orbit[x][54] - 1 * orbit[x][40] - 1 * orbit[x][39] - 1 * orbit[x][32] - 2 * orbit[x][31]) / 3;
      orbit[x][21] = (f_21 - 3 * orbit[x][55] - 3 * orbit[x][50] - 2 * orbit[x][42] - 2 * orbit[x][38] - 2 * orbit[x][33]);
      orbit[x][20] = (f_20 - 2 * orbit[x][54] - 2 * orbit[x][49] - 1 * orbit[x][40] - 1 * orbit[x][37] - 1 * orbit[x][32]);
      orbit[x][19] = (f_19 - 4 * orbit[x][54] - 4 * orbit[x][49] - 1 * orbit[x][40] - 2 * orbit[x][39] - 1 * orbit[x][37] - 2 * orbit[x][35] - 2 * orbit[x][31]);
      orbit[x][18] = (f_18 - 1 * orbit[x][59] - 1 * orbit[x][51] - 2 * orbit[x][46] - 2 * orbit[x][45] - 2 * orbit[x][36] - 2 * orbit[x][27] - 1 * orbit[x][24]) / 2;
      orbit[x][17] = (f_17 - 1 * orbit[x][60] - 1 * orbit[x][53] - 1 * orbit[x][51] - 1 * orbit[x][48] - 1 * orbit[x][37] - 2 * orbit[x][34] - 2 * orbit[x][30]) / 2;
      orbit[x][16] = (f_16 - 1 * orbit[x][59] - 2 * orbit[x][52] - 1 * orbit[x][51] - 2 * orbit[x][46] - 2 * orbit[x][36] - 2 * orbit[x][34] - 1 * orbit[x][29]);
      orbit[x][15] = (f_15 - 1 * orbit[x][59] - 2 * orbit[x][52] - 1 * orbit[x][51] - 2 * orbit[x][45] - 2 * orbit[x][35] - 2 * orbit[x][34] - 2 * orbit[x][27]);

    }
    endTime = clock();
    printf("%.2f sec\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);

    endTime_all = endTime;
    printf("total: %.2f sec\n", (double)(endTime_all - startTime_all) / CLOCKS_PER_SEC);

  } // end count_orbits

  ///
  /// @brief Format the results in an output string
  ///
  std::string format() const
  {
    std::string buffer;
    for (int i = 0; i < this->_n; i++)
    {
      for (int j = 0; j < 73; j++)
      {
        if (j != 0)
        {
          buffer += " ";
        }
        auto orbit_value = this->_orbits[i][j];
        assert(orbit_value >= 0);
        buffer += std::to_string(orbit_value);
      }
      buffer += "\n";
    }
    return buffer;
  } // end count_orbits

}; // end class ORCA


///
/// @brief Intialize an ORCA algorithm from an input graph file
///
class OrcaParser
{
  int _n;
  int _m;
  int _dmax;
  std::vector<key_pair> _edges;
  std::vector<int> _deg;

  ///
  /// @brief Check read values validity
  ///
  void check_validity(int a, int b, int n)
  {
    if (!(0 <= a && a < n) || !(0 <= b && b < n))
    {
      throw std::invalid_argument("Node ids should be between 0 and n-1.");
    }

    if (a == b)
    {
      throw std::invalid_argument("Self loops (edge from x to x) are not allowed.");
    }
  }

  ///
  /// @brief Check for duplicated undirected edges.
  ///
  void check_duplicated_undirected_edges(int m, const std::vector<key_pair> &edges)
  {
    auto first = edges.cbegin();
    auto last = first;
    std::advance(last, m);
    if ( static_cast<int>(std::set<key_pair>(first, last).size()) != m)
    {
      throw std::invalid_argument("Input file contains duplicate undirected edges.");
    }
  }

  ///
  /// @brief Read the number of nodes and edges from file
  ///
  std::istream& read_graph_properties_from(std::ifstream& stream)
  {
    stream >> this->_n >> this->_m;
    return stream;
  }

  ///
  /// @brief Fills the data structures reading from file
  ///
  std::istream& populate_data_from(std::ifstream &stream)
  {
    int a;
    int b;
    int i = 0;
    while(stream >> a >> b)
    {
      check_validity(a, b, this->_n);
      this->_deg[a]++;
      this->_deg[b]++;

      assert( i < this->_edges.size());
      this->_edges.at(i) = key_pair(a, b);
      i=i+1;
    }
    return stream;
  }

  ///
  /// @brief Compute maximum degree of the graph
  ///
  void compute_maximum_degree()
  {
    this-> _dmax = *std::max_element(this->_deg.begin(), this->_deg.end());
  }

public:

  auto n(){return this->_n;}
  auto m(){return this->_m;}
  // Transfer the data, invalidates the parser
  auto take_edges() {return std::move(this->_edges);}
  // Transfer the data, invalidates the parser
  auto take_degrees() {return std::move(this->_deg);}
  // Parse graph from input file
  void parse(const std::string &input_file)
  {
    std::ifstream myfile (input_file);
    if (myfile.is_open())
    {
      read_graph_properties_from(myfile);

      std::cout << "Number of nodes:" << this->_n << std::endl;
      std::cout << "Number of edges:" << this->_m << std::endl;

      this->_edges.resize(this->_m, key_pair(0,0));
      this->_deg.resize(this->_n);

      populate_data_from(myfile);
      compute_maximum_degree();
      check_duplicated_undirected_edges(this->_m, this->_edges);

      std::cout << "Max degree:" << this->_dmax << std::endl;

      myfile.close();

    } else {
      std::cout << "Unable to open file";
    }
  }
  ///
  /// @brief Stream operator can access to private data - not perfect but okay-ish
  ///
  friend std::ostream& operator <<(std::ostream& stream, OrcaParser const& p);
}; // end OrcaParser

///
/// @brief Stream operator
///
std::ostream& operator <<(std::ostream& stream, OrcaParser const& p)
{
  stream << "Number of nodes: " << p._n << std::endl;
  stream << "Number of edges: " << p._m << std::endl;
  stream << "Maximum degree: " << p._dmax << std::endl;
  return stream;
}

///
/// @brief Return a name for the output file based on the input file name
///
/// @note If the input file is "path/to/input.in", the output file should be "path/to/input.out"
///
/// @warning does not work if the filename begins with "."
std::string generate_output_filename_from(const std::string &input_file)
{
  auto rawname = input_file.substr(0, input_file.find_last_of("."));
  return rawname += "_gdvs.out";
}

///
/// @brief Calculate the Graphlet Degree Vector (GDV) for the given graph with the given implementation policy
///
template<typename Policy=adjacent_policy::default_impl>
std::string graphlet_degree_vector_analysis(OrcaParser &parser, const std::string& output_file)
{

  std::cout << "Initializing ORCA algorithm from parser" << std::endl;

  auto computer = ORCA<Policy>
  (
    parser.n(),
    parser.m(),
    std::move(parser.take_edges()),
    std::move(parser.take_degrees())
  );
  // parser has now been emptied from its content

  std::cout << "Counting orbits of graphlets on 5 nodes." << std::endl;
  computer.count_orbits();

  std::cout << "Writing results in " << output_file << std::endl;
  computer.write_results_to(output_file);

  return output_file;
}

///
/// @brief Read data from file, counts orbits, write in output_file
///
std::string read_count_write_orca(const std::string& input_file, const std::string& output_file)
{
  std::cout << "Initialize ORCA from file " << input_file << std::endl;

  OrcaParser parser;
  parser.parse(input_file);
  std::cout << parser << std::endl;

  if( adjacent_policy::compressed::should_use_based_on(parser.n())){
    std::cout << "Using compressed matrix implementation" << std::endl;
    graphlet_degree_vector_analysis<adjacent_policy::compressed>(parser, output_file);
  }else{
    std::cout << "Using default matrix implementation" << std::endl;
    graphlet_degree_vector_analysis<>(parser, output_file);
  }
  return output_file;
}

///
/// @brief Read data from file, counts orbits, write in output_file
///
/// @note output_file defaults to input_file raw name (without extension) + gdvs.out
///
std::string read_count_write_orca(const std::string& input_file)
{
  auto output_file = generate_output_filename_from(input_file);
  read_count_write_orca(input_file, output_file);
  return output_file;
}
#endif
