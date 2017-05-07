#include <bits/stdc++.h>
#include <sys/time.h>
using namespace std;

#define rep(i,n) for(long long i = 0; i < (long long)(n); i++)
#define repi(i,a,b) for(long long i = (long long)(a); i < (long long)(b); i++)
#define pb push_back
#define all(x) (x).begin(), (x).end()
#define fi first
#define se second
#define mt make_tuple
#define mp make_pair
template<class T1, class T2> bool chmin(T1 &a, T2 b) { return b < a && (a = b, true); }
template<class T1, class T2> bool chmax(T1 &a, T2 b) { return a < b && (a = b, true); }
#define exists find_if
#define forall all_of

using ll = long long; using vll = vector<ll>; using vvll = vector<vll>; using P = pair<ll, ll>;
using ld = double;  using vld = vector<ld>; 
using vi = vector<int>; using vvi = vector<vi>; vll conv(vi& v) { vll r(v.size()); rep(i, v.size()) r[i] = v[i]; return r; }
using Pos = complex<double>;

template <typename T, typename U> ostream &operator<<(ostream &o, const pair<T, U> &v) {  o << "(" << v.first << ", " << v.second << ")"; return o; }
template<size_t...> struct seq{}; template<size_t N, size_t... Is> struct gen_seq : gen_seq<N-1, N-1, Is...>{}; template<size_t... Is> struct gen_seq<0, Is...> : seq<Is...>{};
template<class Ch, class Tr, class Tuple, size_t... Is>
void print_tuple(basic_ostream<Ch,Tr>& os, Tuple const& t, seq<Is...>){ using s = int[]; (void)s{0, (void(os << (Is == 0? "" : ", ") << get<Is>(t)), 0)...}; }
template<class Ch, class Tr, class... Args> 
auto operator<<(basic_ostream<Ch, Tr>& os, tuple<Args...> const& t) -> basic_ostream<Ch, Tr>& { os << "("; print_tuple(os, t, gen_seq<sizeof...(Args)>()); return os << ")"; }
ostream &operator<<(ostream &o, const vvll &v) { rep(i, v.size()) { rep(j, v[i].size()) o << v[i][j] << " "; o << endl; } return o; }
template <typename T> ostream &operator<<(ostream &o, const vector<T> &v) { o << '['; rep(i, v.size()) o << v[i] << (i != v.size()-1 ? ", " : ""); o << "]";  return o; }
template <typename T>  ostream &operator<<(ostream &o, const set<T> &m) { o << '['; for (auto it = m.begin(); it != m.end(); it++) o << *it << (next(it) != m.end() ? ", " : ""); o << "]";  return o; }
template <typename T, typename U>  ostream &operator<<(ostream &o, const map<T, U> &m) { o << '['; for (auto it = m.begin(); it != m.end(); it++) o << *it << (next(it) != m.end() ? ", " : ""); o << "]";  return o; }
template <typename T, typename U, typename V>  ostream &operator<<(ostream &o, const unordered_map<T, U, V> &m) { o << '['; for (auto it = m.begin(); it != m.end(); it++) o << *it; o << "]";  return o; }
vector<int> range(const int x, const int y) { vector<int> v(y - x + 1); iota(v.begin(), v.end(), x); return v; }
template <typename T> istream& operator>>(istream& i, vector<T>& o) { rep(j, o.size()) i >> o[j]; return i;}
string bits_to_string(ll input, ll n=64) { string s; rep(i, n) s += '0' + !!(input & (1ll << i)); reverse(all(s)); return s; }

template <typename T> unordered_map<T, ll> counter(vector<T> vec){unordered_map<T, ll> ret; for (auto&& x : vec) ret[x]++; return ret;};
string substr(string s, P x) {return s.substr(x.fi, x.se - x.fi); }
struct ci : public iterator<forward_iterator_tag, ll> { ll n; ci(const ll n) : n(n) { } bool operator==(const ci& x) { return n == x.n; } bool operator!=(const ci& x) { return !(*this == x); } ci &operator++() { n++; return *this; } ll operator*() const { return n; } };

size_t random_seed; namespace std { using argument_type = P; template<> struct hash<argument_type> { size_t operator()(argument_type const& x) const { size_t seed = random_seed; seed ^= hash<ll>{}(x.fi); seed ^= (hash<ll>{}(x.se) << 1); return seed; } }; }; // hash for various class
namespace myhash{ const int Bsizes[]={3,9,13,17,21,25,29,33,37,41,45,49,53,57,61,65,69,73,77,81}; const int xor_nums[]={0x100007d1,0x5ff049c9,0x14560859,0x07087fef,0x3e277d49,0x4dba1f17,0x709c5988,0x05904258,0x1aa71872,0x238819b3,0x7b002bb7,0x1cf91302,0x0012290a,0x1083576b,0x76473e49,0x3d86295b,0x20536814,0x08634f4d,0x115405e8,0x0e6359f2}; const int hash_key=xor_nums[rand()%20]; const int mod_key=xor_nums[rand()%20]; template <typename T> struct myhash{ std::size_t operator()(const T& val) const { return (hash<T>{}(val)%mod_key)^hash_key; } }; };
template <typename T> class uset:public std::unordered_set<T,myhash::myhash<T>> { using SET=std::unordered_set<T,myhash::myhash<T>>; public: uset():SET(){SET::rehash(myhash::Bsizes[rand()%20]);} };
template <typename T,typename U> class umap:public std::unordered_map<T,U,myhash::myhash<T>> { public: using MAP=std::unordered_map<T,U,myhash::myhash<T>>; umap():MAP(){MAP::rehash(myhash::Bsizes[rand()%20]);} };    

struct timeval start; double sec() { struct timeval tv; gettimeofday(&tv, NULL); return (tv.tv_sec - start.tv_sec) + (tv.tv_usec - start.tv_usec) * 1e-6; }
struct init_{init_(){ gettimeofday(&start, NULL); ios::sync_with_stdio(false); cin.tie(0); srand((unsigned int)time(NULL)); random_seed = RAND_MAX / 2 + rand() / 2; }} init__;

static const double EPS = 1e-14;
static const long long INF = 1e18;
static const long long mo = 1e9+7;
#define ldout fixed << setprecision(40) 

// 正規分布の累積確率
// A.m.マーリの方法
double qnorm(double u)
{
    static double a[9] = {   1.24818987e-4, -1.075204047e-3, 5.198775019e-3,
        -0.019198292004, 0.059054035642,-0.151968751364,
        0.319152932694,-0.5319230073,   0.797884560593};
    static double b[15] = { -4.5255659e-5,   1.5252929e-4,  -1.9538132e-5,
        -6.76904986e-4,  1.390604284e-3,-7.9462082e-4,
        -2.034254874e-3, 6.549791214e-3,-0.010557625006,
        0.011630447319,-9.279453341e-3, 5.353579108e-3,
        -2.141268741e-3, 5.35310549e-4,  0.999936657524};
    double w, y, z;
    int i;

    if(u == 0.) return 0.5;
    y = u / 2.;
    if(y < -6.) return 0.;
    if(y > 6.)      return 1.;
    if(y < 0.)      y = - y;
    if(y < 1.) {
        w = y * y;
        z = a[0];
        for(i = 1; i < 9; i++)      z = z * w + a[i];
        z *= (y * 2.);
    } else {
        y -= 2.;
        z = b[0];
        for(i = 1; i < 15; i++) z = z * y + b[i];
    }

    if(u < 0.)  return (1. - z) / 2.;
    return (1. + z) / 2.;
}

// 正規分布のパーセント点
// 戸田の近似式
double pnorm(double qn)
{
    static double b[11] = {  1.570796288,     0.03706987906,  -0.8364353589e-3,
        -0.2250947176e-3, 0.6841218299e-5, 0.5824238515e-5,
        -0.104527497e-5,  0.8360937017e-7,-0.3231081277e-8,
        0.3657763036e-10,0.6936233982e-12};
    double w1, w3;
    int i;

    if(qn < 0. || 1. < qn) {
        fprintf(stderr, "Error : qn <= 0 or qn >= 1  in pnorm()!\n");
        return 0.;
    }
    if(qn == 0.5)   return 0.;

    w1 = qn;
    if(qn > 0.5)    w1 = 1. - w1;
    w3 = -log(4. * w1 * (1. - w1));
    w1 = b[0];
    for(i = 1; i < 11; i++) w1 += (b[i] * pow(w3, (double)i));
    if(qn > 0.5)    return sqrt(w1 * w3);
    return -sqrt(w1 * w3);
}// 自由度nのカイ2乗の上側確率
double qchi(double x2, int n)
{
    double w, pw, x, qc;
    int i, i1;

    if(n < 1) {
        fprintf(stderr,"Error : 自由度 < 1 in qchi()!\n");
        return 0.;
    }
    if(x2 <= 0.)    return 1.;
    if(x2 > 400.)   return 0.;
    if(n > 10) {
        w = 2. / (9. * (double)n);
        pw = pow(x2 / (double)n, 1. / 3.);
        return 1. - qnorm((pw - 1. + w) / sqrt(w));
    }
    w = exp(-x2 / 2.);
    if(n == 2)  return w;
    x = sqrt(x2);
    if(n == 1)  return 2. * (1. - qnorm(x));
    if((n % 2) == 0) {
        i1 = 2;
        qc = w;
    } else {
        i1 = 1;
        qc = 2. * (1. - qnorm(x));
        w *= (0.797884560750774 / x);
    }
    for(i = i1; i <= n - 2; i += 2) {
        w *= (x2 / (double)i);
        qc += w;
    }
    return qc;
}

// 自由度nのカイ2乗のパーセント点
double pchi(double qc, int n)
{
    double c1, c2, gam, x, w, wx;
    int i, j;

    if(qc <= 0. || qc >= 1. || n < 1) {
        fprintf(stderr,"Error : Illigal parameter in pchi()!\n");
        return 0.;
    }
    if(n == 1) {
        w = pnorm(qc / 2.);
        return w * w;
    }
    if(n == 2)  return -2. * log(qc);

    x = -pnorm(qc);
    if(n > 10) {
        w = x * x;
        wx = sqrt(2. * (double)n);
        c1 = (w - 7.) * x / 9. / wx;
        c2 = ((3. * w + 7.) * w - 16.) * 2. / 405. / (double)n;
        wx = (double)n + wx * x + 2. * (w - 1.) / 3. + c1 - c2;
        if(wx < 0.) return 0.;
        return wx;
    }

    w = 2. / 9. / (double)n;
    w = 1. - w + x * sqrt(w);
    wx = (double)n * w * w * w;
    if((n % 2) == 0)    gam = 1.;
    else                gam = 1.772453850905516;
    j = (n + 1) / 2 - 1;
    w = (double)n / 2.;
    for(i = 1; i <= j; i++) gam *= (w - (double)i);
    x = wx / 2.;
    c1 = pow(x, w - 1.);
    c2 = exp(-x) / 2.;
    return wx + (qchi(wx, n) - qc) * gam / c1 / c2;
}

double mean(vector<double> a) {
    double ret = 0; 
    rep(i, a.size()) {
        ret += a[i];
    }
    return ret / a.size();
}
double stdev(vector<double> a) {
    double ret = 0;
    double m = mean(a);
    rep(i, a.size()) {
        ret += (a[i] - m) * (a[i] - m);
    }
    return sqrt(ret / (a.size() - 1));
}

int main(void) {
    vector<vector<double>> a(3, vector<double>(3));
    cin >> a[0][0] >> a[0][1];
    cin >> a[1][0] >> a[1][1];

    rep(i, 2) 
        a[i][2] = a[i][0] + a[i][1];
    rep(i, 2) 
        a[2][i] = a[0][i] + a[1][i];
    a[2][2] = a[0][2] + a[1][2];

    double chisq = 0;
    rep(i, 2) rep(j, 2) {
        double e = a[i][2] * a[2][j] / a[2][2];
        double o_minus_e = a[i][j] - e;
        a[i][j] = o_minus_e * o_minus_e / e;
        chisq += a[i][j];
    }

    if (chisq > pchi(0.05, 1)) {
        cout << "YES" << endl;
    } else {
        cout << "NO" << endl;
    }

    return 0;
}
