#ifndef ORCA_H
#define ORCA_H

#include <string>
#include <algorithm> // std::min, std::max, std::binary_search
#include <unordered_map> // std::unordered_map
#include <fstream> // ofstream
#include <iostream> // std::cout
#include <exception> // std::invalid_argument
#include <functional> // std::function

///
/// @brief Orbit Counting Algorithm
///
/// @note Encapsulate Orca computation logic
///
template<class Strategy>
class ORCA
{
  ///
  /// @brief Nested class for a hashable key for pairs of nodes
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
  /// @brief Nested class for a hashable key for triplet of nodes
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
    key_triple(int a0, int b0, int c0): a(a0), b(b0), c(c0)
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
  namespace adjacent_implementation
  {
    ///
    /// @brief adj is a vector of vector
    ///
    class default
    {
      // adj[x] - adjacency list of node x
      std::vector<std::vector<int>> _adj;
      // Reference on edges
      const std::vector<int>& _deg;
      ///
      /// @brief Checks if an element equivalent to value appears within the range x, y
      ///
      bool operator()(int x, int y) const
      {
        auto first = this->_adj.at(x).cbegin();
        auto last = std::advance(this->_adj.cbegin(), this->_deg.at(x));
        return std::binary_search(first, last, y);
      }
      ///
      /// @brief Build implementation
      ///
      auto build_implementation(int n, int m)
      {
        this->_adj.reserve(n);
        for (int i = 0; i < n; i++)
        {
          _adj[i].reserve(_deg.at(i));
        }
      }
      ///
      /// @brief Policy constructor
      ///
      default(int n, int m, const std::vector<key_pair> &edges, const std::vector<int> & deg):
      _deg(deg),
      _adj(build_implementation(n, m)){}
    };
    ///
    /// @brief adj is one long vector
    ///
    class compressed
    {
      int _n;
      // compressed adjacency matrix
      std::vector<int> _adj;
      // chunk size
      static constexpr int chunk = 8 * sizeof(int);

      ///
      /// @brief Checks if an element equivalent to value appears within the range x, y
      ///
      /// @note Strategy is set up adjacency matrix if it's smaller than 100MB
      ///
      bool operator()(int x, int y) const
      {
        return this->_adj[(x * this->_n + y) / chunk] & (1 << ((x * this->_n + y) % chunk));
      }
      ///
      /// @ brief Build the adjacent matrix
      ///
      auto build_implementation(int n, int m, const std::vector<key_pair> &edges) const
      {
        int size = (n*n) / chunk + 1;
        // Initialize vector
        std::vector<int> matrix.reserve(size);
        // Fill data
        for (int i = 0; i < m; i++)
        {
          int a = edges[i].a;
          int b = edges[i].b;
          _adj[(a * n + b) / chunk] |= (1 << ((a * n + b) % chunk));
          _adj[(b * n + a) / chunk] |= (1 << ((b * n + a) % chunk));
        }
      }
      ///
      /// @brief Policy constructor
      ///
      compressed(int n, int m, const std::vector<key_pair> &edges, const std::vector<int> & deg):
      _n(n),
      _adj_matrix(build_implementation(n, m, edges)){}
    };
  }

  // Type aliases for readibility
  using common2_type = std::unordered_map<key_pair, int, key_pair::hash>;
  using common3_type = std::unordered_map<key_triple, int, key_triple::hash>;

  // ORCA alorithm state (member data)
  common2_type common2;
  // stores the number of nodes that are adjacent to some triplets of nodes
  common3_type common3;

  // number of nodes
  int n;
  // number of edges
  int m;
  // degrees of individual nodes
  std::vector<int> deg;
  // list of edge
  std::vector<key_pair> edges;
  // inc[x] - incidence list of node x: (y, edge id)
  std::vector<std::vector<std::pair<int,int>>> inc;
  // orbit[x][o] - how many times does node x participate in orbit o
  std::vector<std::vector<int>> orbit;

  /// Store a different runtime strategy: can be adjacent_list (default) or adjacent_matrix
  std::function< bool( int, int ) > adjacent;

  ///
  /// @brief Return the value stored at the key, else 0
  ///
  /// @note avoid inflation of the map due to lookups of absent elements.
  ///
  template<class key, class T>
  T find_or_zero(key x, const std::unordered_map<key, T, typename key::hash>& map)
  {
    auto it = map.find(x);
    return (it == map.end()) ? it->second : 0;
  }

  ///
  /// @brief Number of nodes that are adjacent to some pair of nodes
  ///
  int common2_get(const key_pair & x) const
  {
    return find_or_zero(x, this->common2);
  }

  ///
  /// @brief Number of nodes that are adjacent to some triplets of nodes
  ///
  int common3_get(const key_triple & x) const
  {
    return find_or_zero(x, this->common3);
  }

  auto chose_strategy(int n, int m, const& std::vector<key_pair> edges)
  {
    // set up adjacency matrix if it's smaller than 100MB
    if ((int64)n * n < 100LL * 1024 * 1024 * 8)
    {
      // Chose strategy
      this->adjacent = this->adjacent_matrix;

    } else {
      this->adjacent = this->adjacent_list;
    }
  }
  ///
  /// @brief Constructor
  ///
  ORCA(int n, int m, const std::vector<key_pair> &edges, const std::vector<int> &deg):
  adjacent(chose_strategy(n)),

  {



    inc = (PII **)malloc(n * sizeof(PII *));
    for (int i = 0; i < n; i++) {
      inc[i] = (PII *)malloc(deg[i] * sizeof(PII));
    }

    int *d = (int *)calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) {
      int a = edges[i].a, b = edges[i].b;
      adj[a][d[a]] = b;
      adj[b][d[b]] = a;
      inc[a][d[a]] = PII(b, i);
      inc[b][d[b]] = PII(a, i);
      d[a]++;
      d[b]++;
    }

    for (int i = 0; i < n; i++) {
      sort(adj[i], adj[i] + deg[i]);
      sort(inc[i], inc[i] + deg[i]);
    }

    // initialize orbit counts
    orbit = (int64 **)malloc(n * sizeof(int64 *));
    for (int i = 0; i < n; i++) {
      orbit[i] = (int64 *)calloc(73, sizeof(int64));
    }
  }

  ///
  /// @brief Count graphlets on max 5 nodes
  ///
  void count_orbits()
  {
    using PAIR = key_pair;
    using TRIPLE = key_triple;

    clock_t startTime, endTime;
    startTime = clock();

    clock_t startTime_all, endTime_all;
    startTime_all = startTime;

    int frac;
    int frac_prev;

    // precompute common nodes
    printf("stage 1 - precomputing common nodes\n");
    frac_prev = -1;
    for (int x = 0; x < n; x++) {
      frac = 100LL * x / n;
      if (frac != frac_prev) {
        printf("%d%%\r", frac);
        frac_prev = frac;
      }
      for (int n1 = 0; n1 < deg[x]; n1++) {
        int a = adj[x][n1];
        for (int n2 = n1 + 1; n2 < deg[x]; n2++) {
          int b = adj[x][n2];
          PAIR ab = PAIR(a, b);
          common2[ab]++;
          for (int n3 = n2 + 1; n3 < deg[x]; n3++) {
            int c = adj[x][n3];
            int st = adjacent(a, b) + adjacent(a, c) + adjacent(b, c);
            if (st < 2) continue;
            TRIPLE abc = TRIPLE(a, b, c);
            common3[abc]++;
          }
        }
      }
    }

    // precompute triangles that span over edges
    std::vector<int> tri(m);

    for (int i = 0; i < m; i++) {
      int x = edges[i].a, y = edges[i].b;
      for (int xi = 0, yi = 0; xi < deg[x] && yi < deg[y];) {
        if (adj[x][xi] == adj[y][yi]) {
          tri[i]++;
          xi++;
          yi++;
        } else if (adj[x][xi] < adj[y][yi]) {
          xi++;
        } else {
          yi++;
        }
      }
    }

    endTime = clock();
    printf("%.2f sec\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
    startTime = endTime;

    // count full graphlets
    printf("stage 2 - counting full graphlets\n");

    std::vector<int64> C5(n);
    std::vector<int> neigh(n);
    std::vector<int> neigh2(n);

    int nn;
    int nn2;

    // int64 *C5 = (int64 *)calloc(n, sizeof(int64));
    // int *neigh = (int *)malloc(n * sizeof(int)), nn;
    // int *neigh2 = (int *)malloc(n * sizeof(int)), nn2;

    frac_prev = -1;

    for (int x = 0; x < n; x++) {
      frac = 100LL * x / n;
      if (frac != frac_prev) {
        printf("%d%%\r", frac);
        frac_prev = frac;
      }
      for (int nx = 0; nx < deg[x]; nx++) {
        int y = adj[x][nx];
        if (y >= x) break;
        nn = 0;
        for (int ny = 0; ny < deg[y]; ny++) {
          int z = adj[y][ny];
          if (z >= y) break;
          if (adjacent(x, z)) {
            neigh[nn++] = z;
          }
        }
        for (int i = 0; i < nn; i++) {
          int z = neigh[i];
          nn2 = 0;
          for (int j = i + 1; j < nn; j++) {
            int zz = neigh[j];
            if (adjacent(z, zz)) {
              neigh2[nn2++] = zz;
            }
          }
          for (int i2 = 0; i2 < nn2; i2++) {
            int zz = neigh2[i2];
            for (int j2 = i2 + 1; j2 < nn2; j2++) {
              int zzz = neigh2[j2];
              if (adjacent(zz, zzz)) {
                C5[x]++;
                C5[y]++;
                C5[z]++;
                C5[zz]++;
                C5[zzz]++;
              }
            }
          }
        }
      }
    }

    endTime = clock();
    printf("%.2f sec\n", (double)(endTime - startTime) / CLOCKS_PER_SEC);
    startTime = endTime;

    std::vector<int> common_x(n);
    std::vector<int> common_x_list(n);
    int ncx = 0;
    std::vector<int> common_a(n);
    std::vector<int> common_a_list(n);
    int nca = 0;

    // set up a system of equations relating orbit counts
    printf("stage 3 - building systems of equations\n");
    frac_prev = -1;

    for (int x = 0; x < n; x++) {
      frac = 100LL * x / n;
      if (frac != frac_prev) {
        printf("%d%%\r", frac);
        frac_prev = frac;
      }

      for (int i = 0; i < ncx; i++) common_x[common_x_list[i]] = 0;
      ncx = 0;

      // smaller graphlets
      orbit[x][0] = deg[x];
      for (int nx1 = 0; nx1 < deg[x]; nx1++) {
        int a = adj[x][nx1];
        for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
          int b = adj[x][nx2];
          if (adjacent(a, b))
          orbit[x][3]++;
          else
          orbit[x][2]++;
        }
        for (int na = 0; na < deg[a]; na++) {
          int b = adj[a][na];
          if (b != x && !adjacent(x, b)) {
            orbit[x][1]++;
            if (common_x[b] == 0) common_x_list[ncx++] = b;
            common_x[b]++;
          }
        }
      }

      int64 f_71 = 0, f_70 = 0, f_67 = 0, f_66 = 0, f_58 = 0, f_57 = 0;                                // 14
      int64 f_69 = 0, f_68 = 0, f_64 = 0, f_61 = 0, f_60 = 0, f_55 = 0, f_48 = 0, f_42 = 0, f_41 = 0;  // 13
      int64 f_65 = 0, f_63 = 0, f_59 = 0, f_54 = 0, f_47 = 0, f_46 = 0, f_40 = 0;                      // 12
      int64 f_62 = 0, f_53 = 0, f_51 = 0, f_50 = 0, f_49 = 0, f_38 = 0, f_37 = 0, f_36 = 0;            // 8
      int64 f_44 = 0, f_33 = 0, f_30 = 0, f_26 = 0;                                                    // 11
      int64 f_52 = 0, f_43 = 0, f_32 = 0, f_29 = 0, f_25 = 0;                                          // 10
      int64 f_56 = 0, f_45 = 0, f_39 = 0, f_31 = 0, f_28 = 0, f_24 = 0;                                // 9
      int64 f_35 = 0, f_34 = 0, f_27 = 0, f_18 = 0, f_16 = 0, f_15 = 0;                                // 4
      int64 f_17 = 0;                                                                                  // 5
      int64 f_22 = 0, f_20 = 0, f_19 = 0;                                                              // 6
      int64 f_23 = 0, f_21 = 0;                                                                        // 7

      for (int nx1 = 0; nx1 < deg[x]; nx1++) {
        int a = inc[x][nx1].first, xa = inc[x][nx1].second;

        for (int i = 0; i < nca; i++) common_a[common_a_list[i]] = 0;
        nca = 0;
        for (int na = 0; na < deg[a]; na++) {
          int b = adj[a][na];
          for (int nb = 0; nb < deg[b]; nb++) {
            int c = adj[b][nb];
            if (c == a || adjacent(a, c)) continue;
            if (common_a[c] == 0) common_a_list[nca++] = c;
            common_a[c]++;
          }
        }

        // x = orbit-14 (tetrahedron)
        for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first, xb = inc[x][nx2].second;
          if (!adjacent(a, b)) continue;
          for (int nx3 = nx2 + 1; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first, xc = inc[x][nx3].second;
            if (!adjacent(a, c) || !adjacent(b, c)) continue;
            orbit[x][14]++;
            f_70 += common3_get(TRIPLE(a, b, c)) - 1;
            // // debug
            // if (nx2 == nx1 + 6)
            //     printf("%d %d %d\n", x, a, b);
            f_71 += (tri[xa] > 2 && tri[xb] > 2) ? (common3_get(TRIPLE(x, a, b)) - 1) : 0;
            f_71 += (tri[xa] > 2 && tri[xc] > 2) ? (common3_get(TRIPLE(x, a, c)) - 1) : 0;
            f_71 += (tri[xb] > 2 && tri[xc] > 2) ? (common3_get(TRIPLE(x, b, c)) - 1) : 0;
            f_67 += tri[xa] - 2 + tri[xb] - 2 + tri[xc] - 2;
            f_66 += common2_get(PAIR(a, b)) - 2;
            f_66 += common2_get(PAIR(a, c)) - 2;
            f_66 += common2_get(PAIR(b, c)) - 2;
            f_58 += deg[x] - 3;
            f_57 += deg[a] - 3 + deg[b] - 3 + deg[c] - 3;
          }
        }

        // x = orbit-13 (diamond)
        for (int nx2 = 0; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first, xb = inc[x][nx2].second;
          if (!adjacent(a, b)) continue;
          for (int nx3 = nx2 + 1; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first, xc = inc[x][nx3].second;
            if (!adjacent(a, c) || adjacent(b, c)) continue;
            orbit[x][13]++;
            f_69 += (tri[xb] > 1 && tri[xc] > 1) ? (common3_get(TRIPLE(x, b, c)) - 1) : 0;
            f_68 += common3_get(TRIPLE(a, b, c)) - 1;
            f_64 += common2_get(PAIR(b, c)) - 2;
            f_61 += tri[xb] - 1 + tri[xc] - 1;
            f_60 += common2_get(PAIR(a, b)) - 1;
            f_60 += common2_get(PAIR(a, c)) - 1;
            f_55 += tri[xa] - 2;
            f_48 += deg[b] - 2 + deg[c] - 2;
            f_42 += deg[x] - 3;
            f_41 += deg[a] - 3;
          }
        }

        // x = orbit-12 (diamond)
        for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first;
          if (!adjacent(a, b)) continue;
          for (int na = 0; na < deg[a]; na++) {
            int c = inc[a][na].first, ac = inc[a][na].second;
            if (c == x || adjacent(x, c) || !adjacent(b, c)) continue;
            orbit[x][12]++;
            f_65 += (tri[ac] > 1) ? common3_get(TRIPLE(a, b, c)) : 0;
            f_63 += common_x[c] - 2;
            f_59 += tri[ac] - 1 + common2_get(PAIR(b, c)) - 1;
            f_54 += common2_get(PAIR(a, b)) - 2;
            f_47 += deg[x] - 2;
            f_46 += deg[c] - 2;
            f_40 += deg[a] - 3 + deg[b] - 3;
          }
        }

        // x = orbit-8 (cycle)
        for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first, xb = inc[x][nx2].second;
          if (adjacent(a, b)) continue;
          for (int na = 0; na < deg[a]; na++) {
            int c = inc[a][na].first, ac = inc[a][na].second;
            if (c == x || adjacent(x, c) || !adjacent(b, c)) continue;
            orbit[x][8]++;
            f_62 += (tri[ac] > 0) ? common3_get(TRIPLE(a, b, c)) : 0;
            f_53 += tri[xa] + tri[xb];
            f_51 += tri[ac] + common2_get(PAIR(c, b));
            f_50 += common_x[c] - 2;
            f_49 += common_a[b] - 2;
            f_38 += deg[x] - 2;
            f_37 += deg[a] - 2 + deg[b] - 2;
            f_36 += deg[c] - 2;
          }
        }

        // x = orbit-11 (paw)
        for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first;
          if (!adjacent(a, b)) continue;
          for (int nx3 = 0; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first, xc = inc[x][nx3].second;
            if (c == a || c == b || adjacent(a, c) || adjacent(b, c)) continue;
            orbit[x][11]++;
            f_44 += tri[xc];
            f_33 += deg[x] - 3;
            f_30 += deg[c] - 1;
            f_26 += deg[a] - 2 + deg[b] - 2;
          }
        }

        // x = orbit-10 (paw)
        for (int nx2 = 0; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first;
          if (!adjacent(a, b)) continue;
          for (int nb = 0; nb < deg[b]; nb++) {
            int c = inc[b][nb].first, bc = inc[b][nb].second;
            if (c == x || c == a || adjacent(a, c) || adjacent(x, c)) continue;
            orbit[x][10]++;
            f_52 += common_a[c] - 1;
            f_43 += tri[bc];
            f_32 += deg[b] - 3;
            f_29 += deg[c] - 1;
            f_25 += deg[a] - 2;
          }
        }

        // x = orbit-9 (paw)
        for (int na1 = 0; na1 < deg[a]; na1++) {
          int b = inc[a][na1].first, ab = inc[a][na1].second;
          if (b == x || adjacent(x, b)) continue;
          for (int na2 = na1 + 1; na2 < deg[a]; na2++) {
            int c = inc[a][na2].first, ac = inc[a][na2].second;
            if (c == x || !adjacent(b, c) || adjacent(x, c)) continue;
            orbit[x][9]++;
            f_56 += (tri[ab] > 1 && tri[ac] > 1) ? common3_get(TRIPLE(a, b, c)) : 0;
            f_45 += common2_get(PAIR(b, c)) - 1;
            f_39 += tri[ab] - 1 + tri[ac] - 1;
            f_31 += deg[a] - 3;
            f_28 += deg[x] - 1;
            f_24 += deg[b] - 2 + deg[c] - 2;
          }
        }

        // x = orbit-4 (path)
        for (int na = 0; na < deg[a]; na++) {
          int b = inc[a][na].first;
          if (b == x || adjacent(x, b)) continue;
          for (int nb = 0; nb < deg[b]; nb++) {
            int c = inc[b][nb].first, bc = inc[b][nb].second;
            if (c == a || adjacent(a, c) || adjacent(x, c)) continue;
            orbit[x][4]++;
            f_35 += common_a[c] - 1;
            f_34 += common_x[c];
            f_27 += tri[bc];
            f_18 += deg[b] - 2;
            f_16 += deg[x] - 1;
            f_15 += deg[c] - 1;
          }
        }

        // x = orbit-5 (path)
        for (int nx2 = 0; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first;
          if (b == a || adjacent(a, b)) continue;
          for (int nb = 0; nb < deg[b]; nb++) {
            int c = inc[b][nb].first;
            if (c == x || adjacent(a, c) || adjacent(x, c)) continue;
            orbit[x][5]++;
            f_17 += deg[a] - 1;
          }
        }

        // x = orbit-6 (claw)
        for (int na1 = 0; na1 < deg[a]; na1++) {
          int b = inc[a][na1].first;
          if (b == x || adjacent(x, b)) continue;
          for (int na2 = na1 + 1; na2 < deg[a]; na2++) {
            int c = inc[a][na2].first;
            if (c == x || adjacent(x, c) || adjacent(b, c)) continue;
            orbit[x][6]++;
            f_22 += deg[a] - 3;
            f_20 += deg[x] - 1;
            f_19 += deg[b] - 1 + deg[c] - 1;
          }
        }

        // x = orbit-7 (claw)
        for (int nx2 = nx1 + 1; nx2 < deg[x]; nx2++) {
          int b = inc[x][nx2].first;
          if (adjacent(a, b)) continue;
          for (int nx3 = nx2 + 1; nx3 < deg[x]; nx3++) {
            int c = inc[x][nx3].first;
            if (adjacent(a, c) || adjacent(b, c)) continue;
            orbit[x][7]++;
            f_23 += deg[x] - 3;
            f_21 += deg[a] - 1 + deg[b] - 1 + deg[c] - 1;
          }
        }
      }

      // solve equations
      orbit[x][72] = C5[x];
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

  }


  int init(string afin, string afout) {
    // open input, output files
    fin.open(afin.c_str(), fstream::in);
    fout.open(afout.c_str(), fstream::out | fstream::binary);

    if (fin.fail()) {
      cerr << "Failed to open file " << afin << endl;
      return 0;
    }

    if (fout.fail()) {
      cerr << "Failed to open file " << afout << endl;
      return 0;
    }

    // read input graph
    fin >> n >> m;
    int d_max = 0;
    edges = (PAIR *)malloc(m * sizeof(PAIR));
    deg = (int *)calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) {
      int a, b;
      fin >> a >> b;
      if (!(0 <= a && a < n) || !(0 <= b && b < n)) {
        cerr << "Node ids should be between 0 and n-1." << endl;
        return 0;
      }
      if (a == b) {
        cerr << "Self loops (edge from x to x) are not allowed." << endl;
        return 0;
      }
      deg[a]++;
      deg[b]++;
      edges[i] = PAIR(a, b);
    }

    for (int i = 0; i < n; i++) {
      d_max = max(d_max, deg[i]);
    }

    printf("nodes: %d\n", n);
    printf("edges: %d\n", m);
    printf("max degree: %d\n", d_max);
    fin.close();

    if ((int)(set<PAIR>(edges, edges + m).size()) != m) {
      cerr << "Input file contains duplicate undirected edges." << endl;
      return 0;
    }

    // set up adjacency, incidence lists
    adj = (int **)malloc(n * sizeof(int *));
    for (int i = 0; i < n; i++) {
      adj[i] = (int *)malloc(deg[i] * sizeof(int));
    }

    inc = (PII **)malloc(n * sizeof(PII *));
    for (int i = 0; i < n; i++) {
      inc[i] = (PII *)malloc(deg[i] * sizeof(PII));
    }

    int *d = (int *)calloc(n, sizeof(int));
    for (int i = 0; i < m; i++) {
      int a = edges[i].a, b = edges[i].b;
      adj[a][d[a]] = b;
      adj[b][d[b]] = a;
      inc[a][d[a]] = PII(b, i);
      inc[b][d[b]] = PII(a, i);
      d[a]++;
      d[b]++;
    }

    for (int i = 0; i < n; i++) {
      sort(adj[i], adj[i] + deg[i]);
      sort(inc[i], inc[i] + deg[i]);
    }

    // initialize orbit counts
    orbit = (int64 **)malloc(n * sizeof(int64 *));
    for (int i = 0; i < n; i++) {
      orbit[i] = (int64 *)calloc(73, sizeof(int64));
    }

    // Free d
    free(d);

    return 1;
  }

  ///
  /// @brief Format the results in an output string
  ///
  std::string format() const
  {
    std::string buffer;
    for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < 73; j++)
      {
        if (j != 0)
        {
          buffer += " ";
        }
        auto orbit_value = orbit[i][j];
        assert(orbit_value >= 0);
        buffer += std::to_string(orbit_value);
      }
    }
    return buffer
  }

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

}; // end class ORCA

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
  auto last = std::advance(edges.cbegin(), m);
  if ( std::set<key_pair>(first, last).size() != m)
  {
    throw std::invalid_argument("Input file contains duplicate undirected edges.");
  }
}
///
/// @brief Intialize an ORCA algorithm from an input graph file
///
ORCA orca_from_file(const std::string &input_file)
{
  std::string line;
  std::ifstream myfile (input_file);
  int n;
  int m;
  if (myfile.is_open())
  {
    // Read graph properties
    myfile >> n >> m;
    int d_max = 0;
    // Initialize vectors to their size
    std::vector<key_pair> edges(m);
    std::vector<int> deg(n);
    // Populate vectors
    for (int i = 0; i < m; i++)
    {
      int a;
      int b;
      myfile >> a >> b;
      check_validity(a, b, n);
      deg[a]++;
      deg[b]++;
      edges[i] = key_pair(a, b);
    }

    for (int i = 0; i < n; i++)
    {
      d_max = std::max(d_max, deg[i]);
    }

    std::cout << "Number of nodes:" << n << std::endl;
    std::cout << "Number of edges:" << m << std::endl;
    std::cout << "Max degree:" << d_max << std::endl

    myfile.close();

  } else {
    std::cout << "Unable to open file";
  }

  check_duplicated_undirected_edges(m, edges);

  // Call to the constructor
  return ORCA(n, m, edes, deg);
}

///
/// @brief Return a name for the output file based on the input file name
///
/// @note If the input file is "path/to/input.in", the output file should be "path/to/input.out"
///
/// @warning does not work if the filename begins with "."
std::string generate_output_filename(const std::string &input_file)
{
  rawname = input_file.substr(0, input_file.find_last_of("."));
  return rawname += "_gdvs.out";
}

///
/// @brief Calculate the Graphlet Degree Vector (GDV) for the given graph.
///
std::string compute_graphlet_degree_vector(const std::string &input_file)
{
  std::cout << "Initialize ORCA from file " << input_file << std::endl;
  auto computer = orca_from_file(input_file);

  std::cout << "Counting orbits of graphlets on 5 nodes." << std::endl;
  computer.count_orbits();

  auto output_file = generate_output_filename(input_file);

  std::cout << "Writing results in " << output_file << std::endl;
  computer.write_results_to(output_file);

  return output_file;
}

#endif
