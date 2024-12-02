#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include <tuple>

using namespace std;

struct LargeInt {
    bool sign;  // true for positive, false for negative
    vector<uint32_t> LargeIntDigits;
};

const uint32_t RADIX = 1000000000; // Base for storing digits (9 decimal digits)

// Convert hex character to value
uint32_t hexCharToValue(char c) {
    if (isdigit(c)) return c - '0';
    if (isxdigit(c)) return toupper(c) - 'A' + 10;
    throw invalid_argument("Invalid hexadecimal character");
}

// Multiply LargeInt by 16
void multiplyByBase16(LargeInt& number) {
    uint64_t carry = 0;
    for (auto& digit : number.LargeIntDigits) {
        uint64_t result = static_cast<uint64_t>(digit) * 16 + carry;
        digit = result % RADIX;
        carry = result / RADIX;
    }
    if (carry > 0) {
        number.LargeIntDigits.push_back(static_cast<uint32_t>(carry));
    }
}

// Add a single digit to LargeInt
void appendDigit(LargeInt& number, uint32_t digit) {
    uint64_t carry = digit;
    for (auto& element : number.LargeIntDigits) {
        uint64_t sum = static_cast<uint64_t>(element) + carry;
        element = sum % RADIX;
        carry = sum / RADIX;
        if (carry == 0) return;
    }
    if (carry > 0) {
        number.LargeIntDigits.push_back(static_cast<uint32_t>(carry));
    }
}

// Convert hex string to LargeInt
LargeInt hexToLargeInt(const string& hex) {
    LargeInt result = { true, {0} };
    for (char hexChar : hex) {
        multiplyByBase16(result);
        appendDigit(result, hexCharToValue(hexChar));
    }
    return result;
}

// Display LargeInt
void displayLargeInt(const LargeInt& number) {
    if (number.LargeIntDigits.empty()) {
        cout << "0\n";
        return;
    }
    if (!number.sign) cout << "-";

    // Print the most significant digit without leading zeros
    cout << number.LargeIntDigits.back();

    // Print the remaining digits with leading zeros
    for (int i = number.LargeIntDigits.size() - 2; i >= 0; --i) {
        cout << setw(9) << setfill('0') << number.LargeIntDigits[i];
    }
    cout << endl;
}


// Compare two LargeInts
int compareLargeInts(const LargeInt& a, const LargeInt& b) {
    // Check signs first
    if (a.sign != b.sign) {
        return a.sign ? 1 : -1;
    }

    // If signs are the same, compare magnitudes
    if (a.LargeIntDigits.size() != b.LargeIntDigits.size()) {
        return (a.LargeIntDigits.size() > b.LargeIntDigits.size()) ? 1 : -1;
    }

    // Compare digits starting from the most significant
    for (size_t i = a.LargeIntDigits.size(); i-- > 0;) {
        if (a.LargeIntDigits[i] != b.LargeIntDigits[i]) {
            return (a.LargeIntDigits[i] > b.LargeIntDigits[i]) ? 1 : -1;
        }
    }

    // If magnitudes are equal, return 0
    return 0;
}

LargeInt subtractLargeInts(const LargeInt& a, const LargeInt& b);

LargeInt addLargeInts(const LargeInt& a, const LargeInt& b) {
    if (a.sign == b.sign) {
        // If signs are the same, perform addition
        LargeInt result = { a.sign, {} };
        uint64_t carry = 0;
        size_t maxSize = max(a.LargeIntDigits.size(), b.LargeIntDigits.size());
        result.LargeIntDigits.resize(maxSize);

        for (size_t i = 0; i < maxSize; ++i) {
            // Add corresponding digits and carry
            uint64_t sum = carry;
            if (i < a.LargeIntDigits.size()) sum += a.LargeIntDigits[i];
            if (i < b.LargeIntDigits.size()) sum += b.LargeIntDigits[i];

            // Store the result and compute the new carry
            result.LargeIntDigits[i] = sum % RADIX;
            carry = sum / RADIX;
        }

        // Append any remaining carry
        if (carry > 0) {
            result.LargeIntDigits.push_back(static_cast<uint32_t>(carry));
        }
        return result;

    }
    else {
        // If signs differ, convert to subtraction
        if (!a.sign) {
            // a is negative, b is positive => b - |a|
            return subtractLargeInts(b, { true, a.LargeIntDigits });
        }
        else {
            // a is positive, b is negative => a - |b|
            return subtractLargeInts(a, { true, b.LargeIntDigits });
        }
    }
}

LargeInt subtractLargeInts(const LargeInt& a, const LargeInt& b) {
    if (!a.sign && !b.sign) {
        // Both numbers are negative: (-a) - (-b) => b - a
        return subtractLargeInts(b, a);
    }

    if (a.sign != b.sign) {
        // Different signs: a - (-b) => a + b
        return addLargeInts(a, { a.sign, b.LargeIntDigits });
    }

    // Both numbers have the same sign; determine which is larger
    if (compareLargeInts(a, b) < 0) {
        // If |a| < |b|, calculate b - a and flip the sign
        LargeInt result = subtractLargeInts(b, a);
        result.sign = !a.sign; // Negate the sign of the result
        return result;
    }

    // Perform standard subtraction
    LargeInt result = a;
    int64_t borrow = 0;

    for (size_t i = 0; i < b.LargeIntDigits.size() || borrow > 0; ++i) {
        // Subtract the digit and borrow if necessary
        int64_t sub = borrow + (i < b.LargeIntDigits.size() ? b.LargeIntDigits[i] : 0);
        if (result.LargeIntDigits[i] < sub) {
            result.LargeIntDigits[i] += RADIX - sub; // Borrow from the next digit
            borrow = 1;
        }
        else {
            result.LargeIntDigits[i] -= sub;
            borrow = 0;
        }
    }

    // Remove any leading zeros in the result
    while (result.LargeIntDigits.size() > 1 && result.LargeIntDigits.back() == 0) {
        result.LargeIntDigits.pop_back();
    }

    return result;
}

LargeInt multiplyLargeInts(const LargeInt& a, const LargeInt& b) {
    // If either number is zero, return zero
    if (a.LargeIntDigits.empty() || b.LargeIntDigits.empty()) {
        return { true, {0} };
    }

    // Initialize the result with sufficient space for digits
    LargeInt result = { true, vector<uint32_t>(a.LargeIntDigits.size() + b.LargeIntDigits.size(), 0) };

    // Set the sign of the result
    result.sign = (a.sign == b.sign);

    // Perform multiplication digit by digit
    for (size_t i = 0; i < a.LargeIntDigits.size(); ++i) {
        uint64_t carry = 0;
        for (size_t j = 0; j < b.LargeIntDigits.size() || carry > 0; ++j) {
            uint64_t product = result.LargeIntDigits[i + j] + carry;

            if (j < b.LargeIntDigits.size()) {
                product += static_cast<uint64_t>(a.LargeIntDigits[i]) * b.LargeIntDigits[j];
            }

            result.LargeIntDigits[i + j] = product % RADIX;
            carry = product / RADIX;
        }
    }

    // Remove any leading zeros in the result
    while (result.LargeIntDigits.size() > 1 && result.LargeIntDigits.back() == 0) {
        result.LargeIntDigits.pop_back();
    }

    return result;
}

pair<LargeInt, LargeInt> divideLargeInts(const LargeInt& dividend, const LargeInt& divisor) {
    // Handle division by zero
    if (divisor.LargeIntDigits.empty() || (divisor.LargeIntDigits.size() == 1 && divisor.LargeIntDigits[0] == 0)) {
        throw invalid_argument("Division by zero");
    }

    // If dividend is smaller than divisor, quotient = 0, remainder = dividend
    if (compareLargeInts(dividend, divisor) < 0) {
        return { {true, {0}}, dividend };
    }

    // Initialize quotient and remainder
    LargeInt quotient = { dividend.sign == divisor.sign, vector<uint32_t>(dividend.LargeIntDigits.size(), 0) };
    LargeInt remainder = { dividend.sign, dividend.LargeIntDigits };

    // Shift divisor to align with dividend's most significant digit
    size_t shift = dividend.LargeIntDigits.size() - divisor.LargeIntDigits.size();
    LargeInt shiftedDivisor = divisor;
    shiftedDivisor.LargeIntDigits.insert(shiftedDivisor.LargeIntDigits.begin(), shift, 0);

    // Perform division
    for (size_t i = shift + 1; i-- > 0;) {
        uint32_t low = 0, high = RADIX - 1, q = 0;

        // Binary search to find the largest q such that shiftedDivisor * q <= remainder
        while (low <= high) {
            uint32_t mid = low + (high - low) / 2;
            LargeInt testProduct = multiplyLargeInts(shiftedDivisor, { true, {mid} });
            if (compareLargeInts(testProduct, remainder) <= 0) {
                q = mid;
                low = mid + 1;
            }
            else {
                high = mid - 1;
            }
        }

        // Set the quotient and update the remainder
        quotient.LargeIntDigits[i] = q;
        remainder = subtractLargeInts(remainder, multiplyLargeInts(shiftedDivisor, { true, {q} }));

        // Shift the divisor to the right for the next digit
        if (i > 0) {
            shiftedDivisor.LargeIntDigits.erase(shiftedDivisor.LargeIntDigits.begin());
        }
    }

    // Remove leading zeros from quotient and remainder
    while (quotient.LargeIntDigits.size() > 1 && quotient.LargeIntDigits.back() == 0) {
        quotient.LargeIntDigits.pop_back();
    }
    while (remainder.LargeIntDigits.size() > 1 && remainder.LargeIntDigits.back() == 0) {
        remainder.LargeIntDigits.pop_back();
    }

    return { quotient, remainder };
}

LargeInt minLargeInts(const LargeInt& a, const LargeInt& b) {
    return (compareLargeInts(a, b) < 0) ? a : b;
}

// Hàm modLargeInts
LargeInt modLargeInts(const LargeInt& a, const LargeInt& b) {
    if (b.LargeIntDigits.empty() || (b.LargeIntDigits.size() == 1 && b.LargeIntDigits[0] == 0)) {
        throw invalid_argument("Modulo by zero");
    }

    LargeInt remainder = divideLargeInts(a, b).second;

    // Nếu phần dư âm, cộng thêm `b`
    if (!remainder.sign) {
        remainder = addLargeInts(remainder, b);
    }

    return remainder;
}

// Extended Euclidean Algorithm to find gcd and coefficients for LargeInt
tuple<LargeInt, LargeInt, LargeInt> extendedGCD(const LargeInt& a, const LargeInt& b) {
    // Base case
    if (compareLargeInts(b, { true, {0} }) == 0) {
        return make_tuple(a, LargeInt{ true, {1} }, LargeInt{ true, {0} });
    }

    tuple<LargeInt, LargeInt, LargeInt> result = extendedGCD(b, modLargeInts(a, b));  // Recursive call
    LargeInt gcd = get<0>(result);
    LargeInt x1 = get<1>(result);
    LargeInt y1 = get<2>(result);

    // Update coefficients
    LargeInt x = y1;
    LargeInt quotient = divideLargeInts(a, b).first;
    LargeInt y = subtractLargeInts(x1, multiplyLargeInts(quotient, y1));  // Assume LargeInt division and subtraction are implemented

    return make_tuple(gcd, x, y);
}

// Modular Inverse for LargeInt
LargeInt modInverse(const LargeInt& e, const LargeInt& phi) {
    tuple<LargeInt, LargeInt, LargeInt> result = extendedGCD(e, phi);
    LargeInt gcd = get<0>(result);
    LargeInt x = get<1>(result);

    if (compareLargeInts(gcd, { true, {1} }) != 0) {
        throw invalid_argument("e and φ(N) are not coprime, no modular inverse exists.");
    }

    LargeInt modInverse = modLargeInts(x, phi);
    modInverse = addLargeInts(modInverse, phi);
    modInverse = modLargeInts(modInverse, phi);

    // Ensure positive modular inverse
    return modInverse;
}

string LargeIntToHex(const LargeInt& LargeInt) {
    if (LargeInt.LargeIntDigits.empty()) {
        return "0"; // Handle empty LargeInt as zero
    }

    string hexString;

    // Traverse digits in reverse (to handle big-endian output)
    for (int i = LargeInt.LargeIntDigits.size() - 1; i >= 0; --i) {
        stringstream ss;
        ss << hex << uppercase << LargeInt.LargeIntDigits[i];

        string hexPart = ss.str();

        // Ensure each part is padded to 9 hex digits (except the most significant one)
        if (i != LargeInt.LargeIntDigits.size() - 1) {
            while (hexPart.length() < 9) {
                hexPart = "0" + hexPart;
            }
        }

        hexString += hexPart;
    }

    // Remove leading zeros in the final hexadecimal string
    size_t firstNonZero = hexString.find_first_not_of('0');
    if (firstNonZero != string::npos) {
        hexString = hexString.substr(firstNonZero);
    }
    else {
        hexString = "0"; // Handle case where the number is zero
    }

    return hexString;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        cerr << "Usage: " << argv[0] << " <input_file> <output_file>" << endl;
        return 1;
    }

    string inputFile = argv[1];
    string outputFile = argv[2];

    // Read input from the file
    ifstream fin(inputFile);
    if (!fin.is_open()) {
        cerr << "Error: Unable to open input file: " << inputFile << endl;
        return 1;
    }

    string pHex, qHex, eHex;
    getline(fin, pHex);
    getline(fin, qHex);
    getline(fin, eHex);
    fin.close();

    LargeInt p = hexToLargeInt(pHex);
    LargeInt q = hexToLargeInt(qHex);
    LargeInt e = hexToLargeInt(eHex);

    if (compareLargeInts(e, minLargeInts(p, q)) >= 0) {
        cerr << "Error: e must be smaller than min(p, q)" << endl;
        ofstream fout(outputFile);
        if (fout.is_open()) {
            fout << "-1" << endl; // Write error result to the output file
            fout.close();
        }
        return 1;
    }

    LargeInt one = { true, {1} };
    LargeInt pMinus1 = subtractLargeInts(p, one);
    LargeInt qMinus1 = subtractLargeInts(q, one);

    // Compute phi(n) = (p - 1)(q - 1)
    LargeInt phi = multiplyLargeInts(pMinus1, qMinus1);

    // Compute the secret key d
    try {
        LargeInt d = modInverse(e, phi);

        // Write the result to the output file
        ofstream fout(outputFile);
        if (!fout.is_open()) {
            cerr << "Error: Unable to open output file: " << outputFile << endl;
            return 1;
        }

        fout << LargeIntToHex(d) << endl; // Write d in hex format
        fout.close();

    }
    catch (const invalid_argument& ex) {
        cerr << "Error: " << ex.what() << endl;

        // Write error result to the output file
        ofstream fout(outputFile);
        if (fout.is_open()) {
            fout << "-1" << endl;
            fout.close();
        }
        return 1;
    }

    return 0;
}