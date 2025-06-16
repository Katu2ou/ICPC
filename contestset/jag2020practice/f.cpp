//#pragma GCC optimize("O3")
#include <bits/stdc++.h>
using namespace std;

#define rep2(i, m, n) for (int i = (m); i < (n); ++i)
#define rep(i, n) rep2(i, 0, n)
#define drep2(i, m, n) for (int i = (m)-1; i >= (n); --i)
#define drep(i, n) drep2(i, n, 0)
#define all(...) std::begin(__VA_ARGS__), std::end(__VA_ARGS__)
#define INF (long long)1001001001001001001
#define inf 1001001000
#define fi first
#define se second

using ll = long long;
using vi = vector<int>;
using vl = vector<long long>;
using pii = pair<int, int>;
using pll = pair<long long, long long>;
using LL = __int128_t;

ll gcd(ll x, ll y) {if (x == 0) return y;	return gcd(y%x, x);} 
ll lcm(ll x, ll y) { __int128_t xx,yy; xx=x; yy=y; __int128_t ans=xx * yy / gcd(x, y); ll ans2=ans; return ans2; }
template<typename T>
T POW(T x, ll n){T ret=1;	while(n>0){		if(n&1) ret=ret*x;		x=x*x;		n>>=1;	}	return ret;}
template<typename T>
T modpow(T a, ll n, T p) {	if(n==0) return (T)1;  if (n == 1) return a % p;  if (n % 2 == 1) return (a * modpow(a, n - 1, p)) % p;  T t = modpow(a, n / 2, p);  return (t * t) % p;}
template<typename T>
T modinv(T a, T m) {	if(m==0)return (T)1;	T b = m, u = 1, v = 0;	while (b) {		T t = a / b;		a -= t * b; swap(a, b);		u -= t * v; swap(u, v);	}	u %= m;	if (u < 0) u += m;	return u;}
ll REM(ll a, ll b){ return (a % b + b) % b;}
ll QUO(ll a, ll b){ return (a - REM(a, b)) / b;}
/*
random_device rd;
mt19937 gen(rd());
uniform_int_distribution<> dist1(1, 100); [1,100]
int random_value = dist1(gen);
*/
/*
auto start = chrono::steady_clock::now(); //時間計測の開始
auto now = std::chrono::steady_clock::now(); //現在時刻と開始時刻の差を測定
double elapsed = std::chrono::duration<double>(now - start).count(); //時間をdouble型で取得
*/
/*
const int MAXCOMB=510000;
std::vector<mint> FAC(MAXCOMB), FINV(MAXCOMB), INV(MAXCOMB);
void COMinit() {FAC[0] = FAC[1] = 1;FINV[0] = FINV[1] = 1;INV[1] = 1;for (int i = 2; i < MAXCOMB; i++) {FAC[i] = FAC[i - 1] * i;INV[i] = mint(0) - INV[mint::mod() % i] * (mint::mod() / i);FINV[i] = FINV[i - 1] * INV[i];}}
mint COM(int n, int k) {if (n < k) return 0;if (n < 0 || k < 0) return 0;return FAC[n] * FINV[k] * FINV[n - k];}
*/

template <typename T> inline bool chmax(T &a, T b) { return ((a < b) ? (a = b, true) : (false));}
template <typename T> inline bool chmin(T &a, T b) { return ((a > b) ? (a = b, true) : (false));}

/////////////////////////////////////////////////////////////////////////////////////////
int SOLVEFIN = 0;

void solve(){
    for (int i = 0;i < m - 1;i++) {
        ll l = -1, r = x[i] - 1, m;
        if (i > 0) {
            while (r - l > 1) {
                m = (l + r) / 2;
                if ((x[i] - m) * (x[i] - m) < a[i]) r = m;
                else l = m;
            }

            if (r == x[i] - 1) {
                sum[r] = a[i] - 1;
            }
            else if (r == x[i] - 2) {
                sum[r] = a[i] - 4;
                sum[r + 1] = a[i] - 1;
            }
            else {
                b[r] += a[i];
                b[x[i]] -= a[i];

                c[r - 1] -= (x[i] - r) * (x[i] - r);

                e[x[i]] = 1;
                e[r] = -((x[i] - r) * 2 - 1);

                g[x[i]] = 2;
                g[r + 1] = -2;
            }
        }
        l = i, r = n;
        while (r - l > 1) {
            m = (l + r) / 2;
            if ((x[i] - m) * (x[i] - m) < a[i]) l = m;
            else r = m;
        }

        if (l == x[i]) {
            sum[l] = a[i];
        }
        else if (l == x[i] + 1) {
            sum[l] = a[i];
            sum[l + 1] = a[i] - 1;
        }
        else {
            b[x[i]] += a[i];
            b[l + 1] -= a[i];

            d[l + 1] -= (x[i] - l) * (x[i] - l);

            f[x[i] + 1] = 1;
            f[l + 1] = -((l - x[i]) * 2 - 1);

            h[x[i] + 2] = 2;
            h[l + 1] = -2;
        }
    }
    for (int i = 1;i < n;i++) {
        h[i] += h[i - 1];
    }
    for (int i = n - 1;i >= 0;i--) {
        g[i] += g[i + 1];
    }

    for (int i = 1;i < n;i++) {
        f[i] += f[i - 1] + h[i];
    }
    for (int i = n - 1;i >= 0;i--) {
        e[i] += e[i + 1] + g[i];
    }

    for (int i = 1;i < n;i++) {
        d[i] += d[i - 1] + f[i];
    }
    for (int i = n - 1;i >= 0;i--) {
        c[i] += c[i + 1] + e[i];
    }

    for (int i = 1;i < n;i++) b[i] += b[i - 1];
    for (int i = 0;i < n;i++) {
        b[i] -= c[i] + d[i];
    }

    for (int i = 0;i < n;i++) {
        cout << b[i] << endl;
    }
}

signed main(){
	cin.tie(0);
	ios::sync_with_stdio(0);
	cout<<fixed<<setprecision(20);
	while(SOLVEFIN == 0) solve();
}