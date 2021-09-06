#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <sstream>
#include <vector>


uint64_t p;

struct bigint {
public:
    std::vector<uint64_t> numbers;
    uint64_t digit;

    bigint() = default;

    bigint(uint64_t n, uint64_t d): digit(d) {
        while (n > 0) {
            numbers.push_back(n % d);
            n /= d;
        }
        DeleteZeroes();
    }

    bigint(std::vector<uint64_t> v, uint64_t d): numbers(v), digit(d) {
        DeleteZeroes();
    }

    void DeleteZeroes() {
        while (!numbers.empty() && numbers.back() == 0) {
            numbers.pop_back();
        }
    }

    bigint& operator += (const bigint other) {
        uint64_t cur = 0;
        numbers.resize(std::max(other.numbers.size(), numbers.size()) + 1);
        for (uint64_t i = 0; i < numbers.size(); ++i) {
            if (i < other.numbers.size()) {
                numbers[i] += other.numbers[i];
            }
            numbers[i] += cur;
            cur = numbers[i] / digit;
            numbers[i] %= digit;
        }
        DeleteZeroes();
        return *this;
    }

    bigint& operator += (const uint64_t other) {
        *this += bigint(other, digit);
        return *this;
    }

    bigint operator * (const uint64_t other) {
        std::vector<uint64_t> answer;
        uint64_t cur = 0;
        uint64_t i = 0;
        while (cur != 0 || i < numbers.size()) {
            uint64_t res = cur;
            if (i < numbers.size()) {
                res += numbers[i] * other;
            }
            cur = res / digit;
            answer.push_back(res % digit);
            ++i;
        }
        return bigint(answer, digit);
    }

    bigint& operator *= (uint64_t other) {
        *this = *this * other;
        return *this;
    }
};

uint64_t bin_pow(uint64_t pow, uint64_t number) {
    if (pow == 0) {
        return 1;
    }
    if (pow % 2 != 0) {
        return (bin_pow(pow - 1, number) * number) % p;
    }
    auto cur = bin_pow(pow / 2, number);
    return (cur * cur) % p;
}

class Polynomial {
    std::vector<uint64_t> monomials;
    std::vector<uint64_t> module;

public:
    Polynomial() = default;


    Polynomial(const uint64_t lambda, const std::vector<uint64_t> m) {
        monomials = {lambda};
        module = m;
        Norm();
        Remainder();
    }

    Polynomial(const std::vector<uint64_t> v, const std::vector<uint64_t> m)
     : monomials(v), module(m) {
        Norm();
        Remainder();
    }

    std::vector<uint64_t> Module() {
        return module;
    }

    auto Monomials() const {
        return monomials;
    }

    void Remainder() {
        int64_t degree = module.size() - 1;
        while (monomials.size() > degree) {
            uint64_t k = monomials[monomials.size() - 1];
            for (int64_t j = 0; j < degree; ++j) {
                monomials[monomials.size() - j - 2] += p - (k * module[degree - j - 1]) % p;
                monomials[monomials.size() - j - 2] %= p;
            }
            monomials.pop_back();
        }
    }

    friend Polynomial operator * (const Polynomial& p1, const Polynomial& p2) {
        if (p1.Monomials().empty() || p2.Monomials().empty()) {
            return Polynomial(0, p1.module);
        }
        std::vector<uint64_t> res(p1.Monomials().size() + p2.Monomials().size() - 1, 0);
        for (uint64_t i = 0; i < p1.Monomials().size(); ++i) {
            for (uint64_t j = 0; j < p2.Monomials().size(); ++j) {
                res[i + j] += (p1.monomials[i] * p2.monomials[j]) % p;
                res[i + j] %= p;
            }
        }
        return Polynomial(res, p1.module);
    }


    void Norm() {
        if (module.empty()) {
            return;
        }
        if (module.back() != 1) {
            uint64_t k = module.back();
            k = bin_pow(p - 2, k);
            for (auto& c : module) {
                c *= k;
                c %= p;
            }
        }
    }
};

Polynomial bin_pow(uint64_t pow, Polynomial& q) {
    if (pow == 0) {
        std::vector<uint64_t> v = {1};
        return Polynomial(v, q.Module());
    }
    if (pow % 2 != 0) {
        return (bin_pow(pow - 1, q) * q);
    }
    auto cur = bin_pow(pow / 2, q);
    return (cur * cur);
}

char number_to_char(int number) {
    if (0 <= number && number <= 9) {
        return char(number + 48);
    }
    if (10 <= number && number <= 35) {
        return char(number + 55);
    }
    if (36 <= number && number <= 61) {
        return char(number + 61);
    }
    if (number == 62) {
        return char(32);
    }
    if (number == 63) {
        return char(46);
    }
    return char(0);
}


Polynomial string_to_poly(std::string s, std::vector<uint64_t> module) {
    std::stringstream ss(s);
    int64_t x;
    std::vector<uint64_t> vec;
    while (ss >> x) {
        if (x < 0) {
            x += p;
        }
        vec.push_back(x);
    }
    return Polynomial(vec, module);
}

uint64_t Pow(uint64_t pow, uint64_t number) {
    if (pow == 0) {
        return 1;
    }
    if (pow % 2 == 1) {
        return (Pow(pow - 1, number) * number);
    }
    auto cur = Pow(pow / 2, number);
    return (cur * cur);
}

int char_to_number(char symbol) {
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

uint64_t bigint_pow(int n) {
    if (n == 0) {
        return 1;
    }
    if (n % 2 == 1) {
        return bigint_pow(n - 1) * p;
    }
    auto cur = bigint_pow(n / 2);
    return cur * cur;
}

int main() {
    std::cin >> p;
    std::string f_str;
    std::cin.ignore();
    std::getline(std::cin, f_str);
    std::stringstream ss(f_str);
    int64_t x;
    std::vector<uint64_t> f;
    while (ss >> x) {
        if (x < 0) {
            x += p;
        }
        f.push_back(x);
    }
    uint64_t k;
    std::cin >> k;
    std::cin.ignore();
    std::string rs, ms;
    auto Fq = bigint_pow(f.size() - 1);
    std::vector<uint64_t> answer;
    while (std::getline(std::cin, rs)) {
        std::getline(std::cin, ms);
        auto r = string_to_poly(rs, f);
        auto m = string_to_poly(ms, f);
        auto s = bin_pow(Fq - 1 - k, r);
        auto chiper = (m * s);
        for (auto c : chiper.Monomials()) {
            answer.push_back(c);
        }
    }
    bigint number(0, 64);
    for (auto it = answer.rbegin(); it != answer.rend(); ++it) {
        number *= p;
        number += *it;
    }
    for (auto c : number.numbers) {
        std::cout << number_to_char(c);
    }
}