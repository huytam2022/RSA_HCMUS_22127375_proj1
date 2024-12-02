#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <random>

using namespace std;

typedef vector<uint32_t> LargeInt;
constexpr uint32_t RADIX = 1000000000; // 10^9

// Convert hex to value
uint32_t hexCharToValue(char c) {
    if (isdigit(c)) return c - '0';
    if (isxdigit(c)) return toupper(c) - 'A' + 10;
    throw invalid_argument("Invalid hexadecimal character");
}

// Multiply by base 16
void multiplyByBase16(LargeInt& number) {
    uint64_t carryOver = 0;
    for (auto& digit : number) {
        uint64_t result = static_cast<uint64_t>(digit) * 16 + carryOver;
        digit = result % RADIX;
        carryOver = result / RADIX;
    }
    if (carryOver) number.push_back(carryOver);
}

// Add digit to BigInt
void appendDigit(LargeInt& number, uint32_t value) {
    uint64_t carry = value;
    for (auto& digit : number) {
        uint64_t sum = static_cast<uint64_t>(digit) + carry;
        digit = sum % RADIX;
        carry = sum / RADIX;
        if (carry == 0) break;
    }
    if (carry) number.push_back(carry);
}

// Convert hex string to BigInt
LargeInt hexToLargeInt(const string& hex) {
    LargeInt result(1, 0);
    for (char c : hex) {
        multiplyByBase16(result);
        appendDigit(result, hexCharToValue(c));
    }
    return result;
}

// Print BigInt
void displayLargeInt(const LargeInt& number) {
    if (number.empty()) {
        cout << "0\n";
        return;
    }
    cout << number.back();
    for (int i = number.size() - 2; i >= 0; --i) {
        cout << setw(9) << setfill('0') << number[i];
    }
    cout << endl;
}

// Compare two LargeInts
int compareLargeInts(const LargeInt& x, const LargeInt& y) {
    if (x.size() != y.size()) return (x.size() > y.size()) ? 1 : -1;
    for (size_t i = x.size(); i-- > 0;) {
        if (x[i] != y[i]) return (x[i] > y[i]) ? 1 : -1;
    }
    return 0;
}

// Add two BigInt
LargeInt addLargeInts(const LargeInt& a, const LargeInt& b) {
    LargeInt sum;
    size_t maxLength = max(a.size(), b.size());
    sum.resize(maxLength);
    uint64_t carryOver = 0;

    for (size_t i = 0; i < maxLength || carryOver; ++i) {
        if (i == sum.size()) sum.push_back(0);
        uint64_t temp = carryOver;
        if (i < a.size()) temp += a[i];
        if (i < b.size()) temp += b[i];
        sum[i] = temp % RADIX;
        carryOver = temp / RADIX;
    }
    return sum;
}

// Multiply two BigInt
LargeInt multiplyLargeInts(const LargeInt& a, const LargeInt& b) {
    LargeInt product(a.size() + b.size(), 0);
    for (size_t i = 0; i < a.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.size() || carry; ++j) {
            uint64_t temp = product[i + j] + carry;
            if (j < b.size()) temp += static_cast<uint64_t>(a[i]) * b[j];
            product[i + j] = temp % RADIX;
            carry = temp / RADIX;
        }
    }
    while (product.size() > 1 && product.back() == 0) product.pop_back();
    return product;
}

// Subtract two BigInts (a - b)
LargeInt subtractLargeInts(const LargeInt& a, const LargeInt& b) {
    LargeInt result = a;
    int64_t borrow = 0; // Variable to store the borrow for subtraction

    // Traverse through each digit and perform the subtraction
    for (size_t i = 0; i < a.size(); ++i) {
        int64_t sub = borrow + (i < b.size() ? b[i] : 0);

        // If the current digit in 'a' is smaller than the value to subtract, borrow
        if (result[i] < sub) {
            result[i] += RADIX - sub;  // Add BASE (RADIX) and subtract
            borrow = 1;                 // Set borrow for the next iteration
        }
        else {
            result[i] -= sub;          // No borrow, just subtract
            borrow = 0;                // Reset borrow
        }
    }

    // Remove any trailing zeros in the result
    while (result.size() > 1 && result.back() == 0) {
        result.pop_back();
    }

    return result;
}

// Devide two BigInt
pair<LargeInt, LargeInt> divideLargeInts(const LargeInt& dividend, const LargeInt& divisor) {
    if (compareLargeInts(divisor, { 0 }) == 0) throw invalid_argument("Division by zero");
    if (compareLargeInts(dividend, divisor) < 0) return { {0}, dividend };

    size_t shift = dividend.size() - divisor.size();
    LargeInt shiftedDivisor = divisor;
    shiftedDivisor.insert(shiftedDivisor.begin(), shift, 0);
    LargeInt quotient(dividend.size(), 0), remainder = dividend;

    for (size_t i = shift + 1; i-- > 0;) {
        uint32_t low = 0, high = RADIX, currentQuotient = 0;
        while (low <= high) {
            uint32_t mid = low + (high - low) / 2;
            LargeInt temp = multiplyLargeInts(shiftedDivisor, { mid });
            if (compareLargeInts(temp, remainder) <= 0) {
                currentQuotient = mid;
                low = mid + 1;
            }
            else {
                high = mid - 1;
            }
        }
        quotient[i] = currentQuotient;
        remainder = subtractLargeInts(remainder, multiplyLargeInts(shiftedDivisor, { currentQuotient }));
        shiftedDivisor.erase(shiftedDivisor.begin());
    }
    while (quotient.size() > 1 && quotient.back() == 0) quotient.pop_back();
    while (remainder.size() > 1 && remainder.back() == 0) remainder.pop_back();
    return { quotient, remainder };
}

// Modular Exponentiation
LargeInt modExp(LargeInt base, LargeInt exp, const LargeInt& mod) {
    LargeInt result = { 1 };
    base = divideLargeInts(base, mod).second;

    while (!exp.empty() && !(exp.size() == 1 && exp[0] == 0)) {
        if (exp[0] & 1) result = divideLargeInts(multiplyLargeInts(result, base), mod).second;
        base = divideLargeInts(multiplyLargeInts(base, base), mod).second;

        uint32_t carry = 0;
        for (size_t i = exp.size(); i-- > 0;) {
            uint64_t temp = carry * static_cast<uint64_t>(RADIX) + exp[i];
            exp[i] = temp / 2;
            carry = temp % 2;
        }
        while (exp.size() > 1 && exp.back() == 0) exp.pop_back();
    }
    return result;
}

// Main function
int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>\n";
        return 1;
    }

    try {
        ifstream fin(argv[1]);
        ofstream fout(argv[2]);
        if (!fin) throw runtime_error("Cannot open input file");
        if (!fout) throw runtime_error("Cannot open output file");

        string hexNumber;
        fin >> hexNumber;
        LargeInt number = hexToLargeInt(hexNumber);
        bool prime = (number.size() == 1 && number[0] == 2) || modExp({ 2 }, subtractLargeInts(number, { 1 }), number) == LargeInt{ 1 };

        fout << (prime ? 1 : 0) << endl;
    }
    catch (const exception& e) {
        cerr << "Error: " << e.what() << endl;
        return 1;
    }
    return 0;
}
