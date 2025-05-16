# Hoffman encoding/decoding of symbols.
#
# References:
#
# "A Method for the Construction of Minimum-Redundancy Codes", David A. Huffman, Proceedings of the I.R.E., 1952, September
#
# https://videolectures.net/videos/mackay_course_04
#
import heapq
import math
from collections import defaultdict, Counter

FREQ_MAC_KAY5_5 = {"a": 25, "b": 25, "c": 20, "d": 15, "e": 15}
FREQ_MAC_KAY5_7 = {"a": 1, "b": 24, "c": 5, "d": 20, "e": 47, "f": 1, "g": 2}
FREQ_MAC_KAY5_6 = {
    "a": 575,
    "b": 128,
    "c": 263,
    "d": 285,
    "e": 913,
    "f": 173,
    "g": 133,
    "h": 313,
    "i": 599,
    "j": 6,
    "k": 84,
    "l": 335,
    "m": 235,
    "n": 596,
    "o": 689,
    "p": 192,
    "q": 8,
    "r": 508,
    "s": 567,
    "t": 706,
    "u": 334,
    "v": 69,
    "w": 119,
    "x": 73,
    "y": 164,
    "z": 7,
    "–": 1928,
}

FREQ_HUFFMAN_PAPER = {
    "a": 20,
    "b": 18,
    "c": 10,
    "d": 10,
    "e": 10,
    "f": 6,
    "g": 6,
    "h": 4,
    "i": 4,
    "j": 4,
    "k": 4,
    "l": 3,
    "m": 1,
}


def canonical_huffman_codes(symbols_with_lengths):
    """Codes from symbol-length pairs. The pairs must be from the canonical tree constructed by HuffmanNode"""

    # Sort by (length, symbol) as per canonical Huffman rules
    sorted_symbols = sorted(symbols_with_lengths, key=lambda x: (x[1], x[0]))
    code = 0
    prev_len = 0
    codebook = {}
    for symbol, length in sorted_symbols:
        code <<= length - prev_len
        codebook[symbol] = f"{code:0{length}b}"
        code += 1
        prev_len = length
    return codebook


class HuffmanNode:
    def __init__(self, char: chr, freq: int):
        self.freq = freq
        self.char = char
        self.left = None
        self.right = None

    @classmethod
    def from_freq(cls, freq_dict):
        if not freq_dict:
            raise ValueError("Frequency dictionary cannot be empty")
        heap = [cls(char, freq) for char, freq in freq_dict.items()]
        heapq.heapify(heap)
        while len(heap) > 1:
            n1 = heapq.heappop(heap)
            n2 = heapq.heappop(heap)
            merged = cls(char=None, freq=n1.freq + n2.freq)
            merged.left = n1
            merged.right = n2
            heapq.heappush(heap, merged)
        return heap[0]

    @classmethod
    def from_codebook(cls, codebook):
        root = cls(None, 0)

        for symbol, code in codebook.items():
            node = root
            for bit in code:
                if bit == "0":
                    if node.left is None:
                        node.left = cls(None, 0)
                    node = node.left
                else:
                    if node.right is None:
                        node.right = cls(None, 0)
                    node = node.right
            node.char = symbol  # Assign symbol at leaf
        return root

    def display_codes(self):
        codes = self.generate_codes()
        print(f"Codebook ({len(codes)})")
        print("Symbol   l code")
        for c in sorted(codes):
            print(f"{c} -> {len(codes[c]):5} {codes[c]}")

    # def __str__(self):
    #     return f"Node({self.char!r}, freq={self.freq})"

    # def __repr__(self):
    #     return f"HuffmanTree(char={self.char!r}, freq={self.freq})"

    # node comparison - needed for heapq
    def __lt__(self, other):
        return self.freq < other.freq

    def generate_codes(self, prefix="", codebook=None):
        if codebook is None:
            codebook = {}
        if self.char is not None:
            codebook[self.char] = prefix
        else:
            self.left.generate_codes(prefix + "0", codebook)
            self.right.generate_codes(prefix + "1", codebook)
        return codebook

    def encode(self, text):
        codes = self.generate_codes()
        encoded = "".join(codes[c] for c in text)
        return encoded

    def decode(self, encoded):
        if not all(bit in "01" for bit in encoded):
            raise ValueError("Encoded string contains non-binary characters")

        result = []
        node = self
        for bit in encoded:
            node = node.left if bit == "0" else node.right
            if node is None:
                raise ValueError("Invalid encoded string")
            if node.char is not None:
                result.append(node.char)
                node = self
        return "".join(result)


def entropy(freq):
    e = 0.0
    N = sum(freq.values())
    for k in freq:
        if freq[k] > 0:
            p = freq[k] / N
            e -= p * math.log2(p)
    return e


def ex1(text):
    """Calculate symbol frequencies, and create Huffman encoding"""
    freq = Counter(text)
    root = HuffmanNode.from_freq(freq)
    encoded = root.encode(text)
    decoded = root.decode(encoded)

    print("Original:", text)
    print("Encoded:", encoded)
    print("Decoded:", decoded)
    assert text == decoded


def ex2(text):
    root = HuffmanNode.from_freq(FREQ_HUFFMAN_PAPER)
    encoded = root.encode(text)
    decoded = root.decode(encoded)
    print(f"Encoded: {encoded}")
    print("Decoded:", decoded)
    assert text == decoded

    codes = root.generate_codes()
    lav = sum(FREQ_HUFFMAN_PAPER[c] * len(codes[c]) for c in FREQ_HUFFMAN_PAPER)
    assert lav == 342
    root.display_codes()
    print("Average code length: ", lav / sum(FREQ_HUFFMAN_PAPER.values()))


def ex3():
    freq = FREQ_MAC_KAY5_5
    root = HuffmanNode.from_freq(freq)
    codes = root.generate_codes()
    lav = sum(freq[c] * len(codes[c]) for c in freq)
    #root.display_codes()
    #print(f"Average code length: {lav/sum(freq.values())} entropy: {entropy(freq)}")
    assert lav == 230

def ex5():
    freq = FREQ_MAC_KAY5_6
    root = HuffmanNode.from_freq(freq)
    codes = root.generate_codes()
    lav = sum(freq[c] * len(codes[c]) for c in freq)
    root.display_codes()
    print(f"Average code length: {lav/sum(freq.values())} entropy: {entropy(freq)}")
    assert lav == 41462

def ex4():
    # Create Huffman codes, transmit only symbols and code lengths to receiveer
    # which decodes them
    freq = FREQ_MAC_KAY5_7
    root = HuffmanNode.from_freq(freq)
    codes = root.generate_codes()
    symbols_with_lengths = [(c, len(codes[c])) for c in codes]
    codes = canonical_huffman_codes(symbols_with_lengths)
    root = HuffmanNode.from_codebook(codes)
    lav = sum(freq[c] * len(codes[c]) for c in freq)
    assert lav == 197
    root.display_codes()
    print(f"Average code length: {lav/sum(freq.values())} entropy: {entropy(freq)}")


if __name__ == "__main__":
    import sys

    if len(sys.argv) > 1 and sys.argv[1] == "test":
        ex1("this is an example of huffman encoding")
        ex1("abbcccdddd")
        ex2("abbcccdddd")
        ex3()
        ex4()
        ex5()
    else:
        while True:
            print("\nInput a text to be Huffman encoded:")
            text = input()
            freq = Counter(text)
            root = HuffmanNode.from_freq(freq)
            root.display_codes()
            encoded = root.encode(text)
            decoded = root.decode(encoded)

            print("Original:", text)
            print("Encoded:", encoded)
            print("Decoded:", decoded)
            assert text == decoded

            codes = root.generate_codes()
            lav = sum(freq[c] * len(codes[c]) for c in freq)
            lav = lav/sum(freq.values())
            print(f"Average code length: {lav}, Entropy:", entropy(freq))
