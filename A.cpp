#include <algorithm>
#include <deque>
#include <iostream>
#include <random>
#include <string>
#include <vector>


using namespace std;

typedef uint64_t ll;

const uint64_t base = 1000000000;
const uint64_t base_digits = 9;
struct bigint : public bernoulli_distribution {
public:
    vector<uint64_t> a;
    uint64_t sign;
    /*<arpa>*/
    uint64_t size(){
        if(a.empty())return 0;
        uint64_t ans=(a.size()-1)*base_digits;
        uint64_t ca=a.back();
        while(ca)
            ans++,ca/=10;
        return ans;
    }
    bigint operator ^(const bigint &v){
        bigint ans=1,g=*this,b=v;
        while(!b.isZero()){
            if(b%2)
                ans*=g;
            g*=g, b/=2;
        }
        return ans;
    }
    /*</arpa>*/
    bigint() :
            sign(1) {
    }

    bigint(uint64_t v) {
        *this = v;
    }

    bigint(const string &s) {
        read(s);
    }


    void operator=(uint64_t v) {
        sign = 1;
        a.clear();
        if (v < 0)
            sign = -1, v = -v;
        for (; v > 0; v = v / base)
            a.push_back(v % base);
    }

    bigint operator+(const bigint &v) const {
        if (sign == v.sign) {
            bigint res = v;

            for (uint64_t i = 0, carry = 0; i < (uint64_t) std::max(a.size(), v.a.size()) || carry; ++i) {
                if (i == (uint64_t) res.a.size())
                    res.a.push_back(0);
                res.a[i] += carry + (i < (uint64_t) a.size() ? a[i] : 0);
                carry = res.a[i] >= base;
                if (carry)
                    res.a[i] -= base;
            }
            return res;
        }
        return *this - (-v);
    }

    bigint operator-(const bigint &v) const {
        if (sign == v.sign) {
            if (abs() >= v.abs()) {
                bigint res = *this;
                for (uint64_t i = 0, carry = 0; i < (uint64_t) v.a.size() || carry; ++i) {
                    res.a[i] -= carry + (i < (uint64_t) v.a.size() ? v.a[i] : 0);
                    carry = res.a[i] < 0;
                    if (carry)
                        res.a[i] += base;
                }
                res.trim();
                return res;
            }
            return -(v - *this);
        }
        return *this + (-v);
    }

    void operator*=(uint64_t v) {
        if (v < 0)
            sign = -sign, v = -v;
        if(v > base){
            *this = *this * (v / base) * base + *this * (v % base);
            return ;
        }
        for (uint64_t i = 0, carry = 0; i < (uint64_t) a.size() || carry; ++i) {
            if (i == (uint64_t) a.size())
                a.push_back(0);
            uint64_t cur = a[i] * (uint64_t) v + carry;
            carry = (uint64_t) (cur / base);
            a[i] = (uint64_t) (cur % base);
            //asm("divl %%ecx" : "=a"(carry), "=d"(a[i]) : "A"(cur), "c"(base));
        }
        trim();
    }

    bigint operator*(uint64_t v) const {
        bigint res = *this;
        res *= v;
        return res;
    }

    friend pair<bigint, bigint> divmod(const bigint &a1, const bigint &b1) {
        uint64_t norm = base / (b1.a.back() + 1);
        bigint g = a1.abs() * norm;
        bigint b = b1.abs() * norm;
        bigint q, r;
        q.a.resize(g.a.size());

        for (uint64_t i = g.a.size() - 1; i >= 0; i--) {
            r *= base;
            r += g.a[i];
            uint64_t s1 = r.a.size() <= b.a.size() ? 0 : r.a[b.a.size()];
            uint64_t s2 = r.a.size() <= b.a.size() - 1 ? 0 : r.a[b.a.size() - 1];
            uint64_t d = ((uint64_t) base * s1 + s2) / b.a.back();
            r -= b * d;
            while (r < 0)
                r += b, --d;
            q.a[i] = d;
        }

        q.sign = a1.sign * b1.sign;
        r.sign = a1.sign;
        q.trim();
        r.trim();
        return make_pair(q, r / norm);
    }

    bigint operator/(const bigint &v) const {
        return divmod(*this, v).first;
    }

    bigint operator%(const bigint &v) const {
        return divmod(*this, v).second;
    }

    void operator/=(uint64_t v) {
        if (v < 0)
            sign = -sign, v = -v;
        for (uint64_t i = (uint64_t) a.size() - 1, rem = 0; i >= 0; --i) {
            uint64_t cur = a[i] + rem * (uint64_t) base;
            a[i] = (uint64_t) (cur / v);
            rem = (uint64_t) (cur % v);
            if (i == 0) {
                break;
            }
        }
        trim();
    }

    bigint operator/(uint64_t v) const {
        bigint res = *this;
        res /= v;
        return res;
    }

    uint64_t operator%(uint64_t v) const {
        uint64_t m = 0;
        for (uint64_t i = a.size() - 1; i >= 0; --i) {
            m = (a[i] + (m * ((uint64_t) base % v) % v)) % v;
            if (i == 0) {
                break;
            }
        }
        return m * sign;
    }

    void operator+=(const bigint &v) {
        *this = *this + v;
    }
    void operator-=(const bigint &v) {
        *this = *this - v;
    }
    void operator*=(const bigint &v) {
        *this = *this * v;
    }
    void operator/=(const bigint &v) {
        *this = *this / v;
    }

    bool operator<(const bigint &v) const {
        if (sign != v.sign)
            return sign < v.sign;
        if (a.size() != v.a.size())
            return a.size() * sign < v.a.size() * v.sign;
        for (int i = a.size() - 1; i >= 0; i--)
            if (a[i] != v.a[i])
                return a[i] * sign < v.a[i] * sign;
        return false;
    }

    bool operator>(const bigint &v) const {
        return v < *this;
    }
    bool operator<=(const bigint &v) const {
        return !(v < *this);
    }
    bool operator>=(const bigint &v) const {
        return !(*this < v);
    }
    bool operator==(const bigint &v) const {
        return !(*this < v) && !(v < *this);
    }
    bool operator!=(const bigint &v) const {
        return *this < v || v < *this;
    }

    void trim() {
        while (!a.empty() && !a.back())
            a.pop_back();
        if (a.empty())
            sign = 1;
    }

    bool isZero() const {
        return a.empty() || (a.size() == 1 && !a[0]);
    }

    bigint operator-() const {
        bigint res = *this;
        res.sign = -sign;
        return res;
    }

    bigint abs() const {
        bigint res = *this;
        res.sign *= res.sign;
        return res;
    }

    uint64_t longValue() const {
        uint64_t res = 0;
        for (uint64_t i = a.size() - 1; i >= 0; i--)
            res = res * base + a[i];
        return res * sign;
    }

    friend bigint gcd(const bigint &g, const bigint &b) {
        return b.isZero() ? g : gcd(b, g % b);
    }
    friend bigint lcm(const bigint &g, const bigint &b) {
        return g / gcd(g, b) * b;

    }

    void read(const string &s) {
        sign = 1;
        a.clear();
        uint64_t pos = 0;
        while (pos < (uint64_t) s.size() && (s[pos] == '-' || s[pos] == '+')) {
            if (s[pos] == '-')
                sign = -sign;
            ++pos;
        }
        for (uint64_t i = s.size() - 1; i >= pos; i -= base_digits) {
            uint64_t x = 0;
            for (uint64_t j = std::max(pos, i - base_digits + 1); j <= i; j++)
                x = x * 10 + s[j] - '0';
            a.push_back(x);
        }
        trim();
    }

    friend istream& operator>>(istream &stream, bigint &v) {
        string s;
        stream >> s;
        v.read(s);
        return stream;
    }

    static vector<uint64_t> convert_base(const vector<uint64_t> &a, uint64_t old_digits, uint64_t new_digits) {
        vector<uint64_t> p(std::max(old_digits, new_digits) + 1);
        p[0] = 1;
        for (uint64_t i = 1; i < (uint64_t) p.size(); i++)
            p[i] = p[i - 1] * 10;
        vector<uint64_t> res;
        uint64_t cur = 0;
        uint64_t cur_digits = 0;
        for (uint64_t i = 0; i < (uint64_t) a.size(); i++) {
            cur += a[i] * p[cur_digits];
            cur_digits += old_digits;
            while (cur_digits >= new_digits) {
                res.push_back(uint64_t(cur % p[new_digits]));
                cur /= p[new_digits];
                cur_digits -= new_digits;
            }
        }
        res.push_back((uint64_t) cur);
        while (!res.empty() && !res.back())
            res.pop_back();
        return res;
    }

    typedef vector<uint64_t> vll;

    static vll karatsubaMultiply(const vll &a, const vll &b) {
        uint64_t n = a.size();
        vll res(n + n);
        if (n <= 32) {
            for (uint64_t i = 0; i < n; i++)
                for (uint64_t j = 0; j < n; j++)
                    res[i + j] += a[i] * b[j];
            return res;
        }

        uint64_t k = n >> 1;
        vll a1(a.begin(), a.begin() + k);
        vll a2(a.begin() + k, a.end());
        vll b1(b.begin(), b.begin() + k);
        vll b2(b.begin() + k, b.end());

        vll a1b1 = karatsubaMultiply(a1, b1);
        vll a2b2 = karatsubaMultiply(a2, b2);

        for (uint64_t i = 0; i < k; i++)
            a2[i] += a1[i];
        for (uint64_t i = 0; i < k; i++)
            b2[i] += b1[i];

        vll r = karatsubaMultiply(a2, b2);
        for (uint64_t i = 0; i < (uint64_t) a1b1.size(); i++)
            r[i] -= a1b1[i];
        for (uint64_t i = 0; i < (uint64_t) a2b2.size(); i++)
            r[i] -= a2b2[i];

        for (uint64_t i = 0; i < (uint64_t) r.size(); i++)
            res[i + k] += r[i];
        for (uint64_t i = 0; i < (uint64_t) a1b1.size(); i++)
            res[i] += a1b1[i];
        for (uint64_t i = 0; i < (uint64_t) a2b2.size(); i++)
            res[i + n] += a2b2[i];
        return res;
    }

    bigint operator*(const bigint &v) const {
        vector<uint64_t> a6 = convert_base(this->a, base_digits, 6);
        vector<uint64_t> b6 = convert_base(v.a, base_digits, 6);
        vll g(a6.begin(), a6.end());
        vll b(b6.begin(), b6.end());
        while (g.size() < b.size())
            g.push_back(0);
        while (b.size() < g.size())
            b.push_back(0);
        while (g.size() & (g.size() - 1))
            g.push_back(0), b.push_back(0);
        vll c = karatsubaMultiply(g, b);
        bigint res;
        res.sign = sign * v.sign;
        for (uint64_t i = 0, carry = 0; i < (uint64_t) c.size(); i++) {
            uint64_t cur = c[i] + carry;
            res.a.push_back((uint64_t) (cur % 1000000));
            carry = (uint64_t) (cur / 1000000);
        }
        res.a = convert_base(res.a, 6, base_digits);
        res.trim();
        return res;
    }
};

uint64_t char_to_number(char symbol) {
    if (symbol >= 48 && symbol <= 57)
        return symbol - 48;
    if (symbol >= 65 && symbol <= 90)
        return symbol - 55;
    if (symbol >= 97 && symbol <= 122)
        return symbol - 61;
    if (symbol == 32)
        return 62;
    if (symbol == 46)
        return 63;
    return 64;
}

uint64_t bin_pow(uint64_t p, uint64_t pow, uint64_t number) {
    if (pow == 0) {
        return 1;
    }
    if (pow % 2 == 1) {
        return (bin_pow(p, pow - 1, number) * number) % p;
    }
    auto cur = bin_pow(p, pow / 2, number);
    return (cur * cur) % p;
}

int main() {
    std::ios_base::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::mt19937 rand;
    uint64_t g, k, p;
    std::cin >> p >> g >> k;
    bigint number(0);
    std::string s;
    std::getline(std::cin, s);
    std::getline(std::cin, s);
    for (auto it = s.rbegin(); it != s.rend(); ++it) {
        number *= 64;
        number += char_to_number(*it);
    }
    while (number > 0) {
        uint64_t digit = number % p;
        uint64_t b = rand() % (p - 1) + 1;
        uint64_t ans2 = (bin_pow(p, b, k) * digit) % p;
        uint64_t ans1 = bin_pow(p, b, g) % p;
        std::cout << ans1 << " " << ans2 << "\n";
        number /= p;
    }
}
