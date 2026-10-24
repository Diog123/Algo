// Description: 
// hash(s[l..r]) = sum_{k=l}^{r} val(s[k]) * p^(k - l)   (mod m)
//
// getHash(l, r) is normalized: it returns exactly what hashing
// s.substr(l, r - l + 1) from scratch would give, so hashes can be compared
// across different positions and different strings.
// basically we divide with inverse and we don't use the multiplication trick.
class hashedString {
public:
	hashedString() : hashedString(std::string{}) {}

	hashedString(std::string s) : str(s) {
		strSize = str.size();
		build();
	}

	std::size_t length() const {
		return strSize;
	}

	int getHash() const { // only the first type
		return H[0].back();
	}

	int getHash(int l, int r) const { // s[l..r] only the first type.
		return getRangeHash(l, r, 0);
	}

	std::vector<int> getHashes() const { // if hash collisions this is safer
		return getHashes(0, static_cast<int>(strSize) - 1);
	}

	std::vector<int> getHashes(int l, int r) const {
		std::vector<int> result(modulos.size());
		for (std::size_t htype = 0; htype < modulos.size(); htype++) {
			result[htype] = getRangeHash(l, r, htype);
		}
		return result;
	}

private:
	std::string str{};
	// It's best to keep the moduli in int, so every product fits in 64 bits.
	// Every modulus must be prime: p^(-1) is computed with Fermat's little theorem.
	std::vector<int> modulos = {1'000'000'007, 1'000'000'009};
	std::vector<int> powers = {31, 33};

	std::vector<std::vector<int>> H;    // H[t][i]    = hash of the prefix s[0..i-1]
	std::vector<std::vector<int>> P;    // P[t][i]    = p^i    (mod m)
	std::vector<std::vector<int>> invP; // invP[t][i] = p^(-i) (mod m)
	std::size_t strSize{};

	static int charValue(char c) {
		return c - 'a' + 1;
	}

	static std::int64_t modPow(std::int64_t base, std::int64_t exp, std::int64_t mod) {
		std::int64_t result = 1 % mod;
		base %= mod;
		while (exp > 0) {
			if (exp & 1) {
				result = result * base % mod;
			}
			base = base * base % mod;
			exp >>= 1;
		}
		return result;
	}

	void build() {
		assert(!modulos.empty() && modulos.size() == powers.size());
		const std::size_t hashCount = modulos.size();

		H.assign(hashCount, std::vector<int>(strSize + 1));
		P.assign(hashCount, std::vector<int>(strSize + 1));
		invP.assign(hashCount, std::vector<int>(strSize + 1));

		for (std::size_t htype = 0; htype < hashCount; htype++) {
			const std::int64_t m = modulos[htype];
			const std::int64_t p = powers[htype];

			const std::int64_t pInv = modPow(p, m - 2, m);
			
			H[htype][0] = 0;
			P[htype][0] = 1;
			invP[htype][0] = 1;
			for (std::size_t i = 1; i <= strSize; i++) {
				P[htype][i] = static_cast<int>(P[htype][i - 1] * p % m);
				invP[htype][i] = static_cast<int>(invP[htype][i - 1] * pInv % m);
				// s[i-1] gets weight p^(i-1), i.e. hash(s) = sum s[k] * p^k
				const std::int64_t term = std::int64_t{charValue(str[i - 1])} * P[htype][i - 1] % m;
				H[htype][i] = static_cast<int>((H[htype][i - 1] + term) % m);
			}
		}
	}

	// Hash of s[l..r] for one hash type, normalized to start at weight p^0.
	int getRangeHash(int l, int r, std::size_t htype) const {
		const std::int64_t m = modulos[htype];
		// H[r+1] - H[l] = hash(s[l..r]) * p^l   (mod m)
		const std::int64_t shifted = (H[htype][r + 1] - H[htype][l] + m) % m;
		// dividing by p^l is the same as multiplying by p^(-l)
		return static_cast<int>(shifted * invP[htype][l] % m);
	}
};
