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
}
// 自由度nのt分布の両側確率
void qtsub(double *q, int n, double w2, double w3, double t2)
{
    int i, j;
    double w;

    j = (n - 2) / 2;
    for(i = 1; i <= j; i++) {
        w = 2. * (double)i - w3;
        *q -= (w2 *= (t2 * w / (1. + w)));
    }
}
double qt(double t, int n)
{
    double t1, t2, w1, w2, wq;

    if(n < 1)
    {
        fprintf(stderr,"Error : n < 1  in qt()!\n");
        return 0.;
    }

    w1 = 0.636619772284456;
    if(t < 0.)  t = - t;
    t1 = t / sqrt((double)n);
    t2 = 1. / (1. + t1 * t1);
    if((n % 2) != 0) {
        wq = 1. - w1 * atan(t1);
        if(n != 1) {
            wq -= (w2 = w1 * t1 * t2);
            if (n != 3) qtsub(&wq, n, w2, 0., t2);
        }
        if(wq > 0.) return wq;
        return 0.;
    }
    wq = 1. - (w2 = t1 * sqrt(t2));
    if(n != 2)      qtsub(&wq, n, w2, 1., t2);
    if(wq > 0.) return wq;
    return 0.;
}


// 自由度nのt分布のパーセント点
double ptsub(double q, int n)
{
    double eps, qe, s, w;

    if(n == 1 && 0.001 <= q && q < 0.01)    eps = 1.e-4;
    else if (n == 2 && q < 0.0001)          eps = 1.e-4;
    else if (n == 1 && q < 0.001)           eps = 1.e-2;
    else                                    eps = 1.e-5;
    s = 10000.;
    w = 0.;
    for(;;) {
        w += s;
        if(s <= eps)    return w;
        if((qe = qt(w,n) - q) == 0.)    return w;
        if(qe < 0.) {
            w -= s;
            s /= 10.;
        }
    }
} 
double pt(double q, int n)
{
    double f, f1, f2, f3, f4, f5, u, u2, w, w0, w1, w2, w3, w4;

    if(q < 1.e-5 || q > 1. || n < 1) {
        fprintf(stderr,"Error : Illigal parameter  in pt()!\n");
        return 0.;
    }

    if(n <= 5)  return ptsub(q, n);

    if(q <= 5.e-3 && n <= 13)   return ptsub(q, n);

    f1 = 4. * (f = (double)n);
    f5 = (f4 = (f3 = (f2 = f * f) * f) * f) * f;
    f2 *= 96.;
    f3 *= 384.;
    f4 *= 92160.;
    f5 *= 368640.;
    u = pnorm(1. - q / 2.);

    w0 = (u2 = u * u) * u;
    w1 = w0 * u2;
    w2 = w1 * u2;
    w3 = w2 * u2;
    w4 = w3 * u2;
    w = ((w0 + u) / f1);
    w += ((5. * w1 + 16. * w0 + 3. * u) / f2);
    w += ((3. * w2 + 19. * w1 + 17. * w0 - 15. * u) / f3);
    w += ((79. * w3 + 776. * w2 + 1482. * w1 - 1920. * w0 - 945. * u) / f4);
    w += ((27. * w4 + 339. * w3 + 930. * w2 - 1782. * w1 - 765. * w0
                + 17955. * u) / f5);
    return u + w;
}

// 自由度n1, n2のF分布の上側確率
double qf(double f, int n1, int n2)
{
    double d, p, w, y, z;
    int j, k, m, mn, n;

    if(f < 0. || n1 < 1 || n2 < 1) {
        fprintf(stderr,"Error : Illegal parameter  in qf()!\n");
        return 0.;
    }

    w = f * (double)n1 / (double)n2;
    z = 1. / (1. + w);
    m = 2 * (n1 / 2) - n1 + 2;
    n = 2 * (n2 / 2) - n2 + 2;
    if((mn = 4 - 2 * (n1 % 2) - (n2 % 2)) == 1) {
        y = 0.318309886142228;
        p = sqrt(w);
        d = y * z / p;
        p = 2. * y * atan(p);
    } else if(mn == 2) {
        p = sqrt(z * w);
        d = p * z / w / 2.;
    } else if(mn == 3) {
        p = sqrt(z);
        d = z * p / 2.;
        p = 1. - p;
    } else {
        d = z * z;
        p = z * w;
    }
    if(n1 <= 2 && n2 <= 2)  return 1. - p;

    y = 2. * w / z;
    if((k = n + 2) <= n2) {
        for( ; k <= n2; k += 2) {
            d *= ((1. + (double)m / (double)(k - 2)) * z);
            if(m != 1)  p = (p + w) * z;
            else        p += (d * y / (double)(k - 1));
        }
    }
    k = m + 2;
    if(k > n1)  return 1. - p;

    y = z * w;
    z = 2. / z;
    j = n2 - 2;
    for( ; k <= n1; k += 2) {
        d *= (y * (double)(k + j) / (double)(k - 2));
        p -= (z * d / (double)(k + j));
    }
    return 1. - p;
}

double pfsub(double x, double y, double z) {
    return (sqrt(z) - y) / x / 2.;
}

// 自由度n1, n2のF分布のパーセント点
double pf(double q, int n1, int n2)
{
    double a, b, c, d, eps, fw, qe, qn, s, u, u2, w1, w2, w3, w4;

    if(q < 0. || q > 1. || n1 < 1 || n2 < 1) {
        fprintf(stderr,"Error : Illegal parameter  in pf()!\n");
        return 0.;
    }

    if(n1 <= 240 || n2 <= 240) {
        eps = 1.e-5;
        if(n2 == 1) eps = 1.e-4;
        s = 1000.;
        fw = 0.;
        for(;;) {
            fw += s;
            if(s <= eps)    return fw;
            if((qe = qf(fw,n1,n2) - q) == 0.)   return fw;
            if(qe < 0.) {
                fw -= s;
                s /= 10.;
            }
        }
    }

    eps = 1.e-6;
    qn = q;
    if(q < 0.5) qn = 1. - q;
    u = pnorm(qn);
    w1 = 2. / (double)n1 / 9.;
    w2 = 2. / (double)n2 / 9.;
    w3 = 1. - w1;
    w4 = 1. - w2;
    u2 = u * u;
    a = w4 * w4 - u2 * w2;
    b = -2. * w3 * w4;
    c = w3 * w3 - u2 * w1;
    d = b * b - 4 * a * c;
    if(d < 0.)  
        fw = pfsub(a, b, 0.);
    else {
        if(fabs(a) > eps)
            fw = pfsub(a, b, d);
        else {
            if(fabs(b) > eps)   return -c / b;
            fw = pfsub(a, b, 0.);
        }
    }
    return fw * fw * fw;
}
// http://www5.airnet.ne.jp/tomy/cpro/sslib11.htm



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
    int m, n; cin >> m >> n;
    vector<double> a(m), b(n); cin >> a >> b;
    
    // (m, a, sa)のほうが分散が多くswap
    double sa = stdev(a), sb = stdev(b);
    if (sa < sb) swap(sa, sb), swap(n, m), swap(a, b);
    double ua = sa * sa, ub = sb * sb;

    double f = ua / ub;
    double t = 0, v = 0;
    if (f < pf(0.05, m-1, n-1)) {
        // 等分散仮定OK
        // 等分散のt検定, 両側
        double u = ((m-1)*ua+(n-1)*ub) / (n + m - 2);
        t = abs(mean(a) - mean(b)) / sqrt(u * (1.0 / m + 1.0 / n));
        v = n + m - 2;
    } else {
        // 等分散仮定not OK
        // ウェルチのt検定, 両側
        t = abs(mean(a) - mean(b)) / sqrt(ua/m+ub/n);
        v = pow(ua/m+ub/n, 2) / (ua*ua/m/m/(m-1)+ub*ub/n/n/(n-1));
    }
    if (t > pt(0.05, v)){
        cout << "YES" << endl;
    } else {
        cout << "NO" << endl;
    }
    return 0;
}
