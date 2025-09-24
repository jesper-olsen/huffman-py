# Huffman encoding/decoding of symbols.
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


def canonical_huffman_codes(symbols_with_lengths):
    """
    Generate canonical Huffman codes from symbol-length pairs.
    
    Args:
        symbols_with_lengths: List of (symbol, code_length) tuples
        
    Returns:
        dict: Mapping of symbols to their binary code strings
    """

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
    def __init__(self, char, freq: int):
        self.freq = freq
        self.char = char
        self.left = None
        self.right = None
        self.codebook = None

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
        heap[0].generate_codes()
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
        root.codebook = codebook
        return root

    def display_codes(self, freq):
        print(f"Codebook ({len(self.codebook)})")
        print("Symbol  cnt length  code")
        for c in sorted(self.codebook):
            print(f"{c} ->   {freq[c]:4}  {len(self.codebook[c]):5}  {self.codebook[c]}")

    def __str__(self):
        if self.char is not None:
            return f"Leaf({self.char!r}, freq={self.freq})"
        return f"Internal(freq={self.freq})"

    def __repr__(self):
        if self.char is not None:
            return f"HuffmanNode(char={self.char!r}, freq={self.freq})"
        return f"HuffmanNode(freq={self.freq}, left={bool(self.left)}, right={bool(self.right)})"

    def pretty_print(self, indent=0):
        prefix = "  " * indent
        if self.char is not None:
            print(f"{prefix}Leaf('{self.char}', freq={self.freq})")
        else:
            print(f"{prefix}Node(freq={self.freq})")
            if self.left:
                self.left.pretty_print(indent + 1)
            if self.right:
                self.right.pretty_print(indent + 1)

    # node comparison - needed for heapq
    def __lt__(self, other):
        return self.freq < other.freq

    def generate_codes(self):
        if self.codebook is not None:
            return self.codebook

        self.codebook = {}
        if self.char is not None:
            # Edge case: tree with single node
            self.codebook[self.char] = "0"
        else:
            self._build_codes("", self.codebook)

    def _build_codes(self, prefix, codebook):
        if self.char is not None:
            codebook[self.char] = prefix
        else:
            if self.left:
                self.left._build_codes(prefix + "0", codebook)
            if self.right:
                self.right._build_codes(prefix + "1", codebook)

    def encode(self, text):
        missing_chars = set(text) - set(self.codebook.keys())
        if missing_chars:
            raise ValueError(f"Characters not in tree: {missing_chars}")
        encoded = "".join(self.codebook[c] for c in text)
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

        if node != self:
            raise ValueError("Incomplete encoded string - ended mid-traversal")

        return "".join(result)


def entropy(freq):
    N = sum(freq.values())
    if N == 0:
        return 0.0
    return -sum((f / N) * math.log2(f / N) for f in freq.values() if f > 0)


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
    freq = FREQ_HUFFMAN_PAPER
    root = HuffmanNode.from_freq(freq)
    encoded = root.encode(text)
    decoded = root.decode(encoded)
    print(f"Encoded: {encoded}")
    print("Decoded:", decoded)
    assert text == decoded

    lav = sum(freq[c] * len(root.codebook[c]) for c in freq)
    assert lav == 342
    root.display_codes(freq)
    print("Average code length: ", lav / sum(freq.values()))


def ex3():
    FREQ_MAC_KAY5_5 = {"a": 25, "b": 25, "c": 20, "d": 15, "e": 15}
    freq = FREQ_MAC_KAY5_5
    root = HuffmanNode.from_freq(freq)
    lav = sum(freq[c] * len(root.codebook[c]) for c in freq)
    # root.display_codes()
    # print(f"Average code length: {lav/sum(freq.values())} entropy: {entropy(freq)}")
    assert lav == 230


def ex5():
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
    freq = FREQ_MAC_KAY5_6
    root = HuffmanNode.from_freq(freq)
    lav = sum(freq[c] * len(root.codebook[c]) for c in freq)
    root.display_codes(freq)
    print(f"Average code length: {lav/sum(freq.values())} entropy: {entropy(freq)}")
    assert lav == 41462


def ex4():
    # Create Huffman codes, transmit only symbols and code lengths to receiveer
    # which decodes them
    FREQ_MAC_KAY5_7 = {"a": 1, "b": 24, "c": 5, "d": 20, "e": 47, "f": 1, "g": 2}
    freq = FREQ_MAC_KAY5_7
    root = HuffmanNode.from_freq(freq)
    symbols_with_lengths = [(c, len(v)) for c,v in root.codebook.items()]
    codes = canonical_huffman_codes(symbols_with_lengths)
    root = HuffmanNode.from_codebook(codes)
    lav = sum(freq[c] * len(codes[c]) for c in freq)
    assert lav == 197
    root.display_codes(freq)
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
            encoded = root.encode(text)
            decoded = root.decode(encoded)

            print("Encoded:", encoded)
            print("Decoded:", decoded)
            assert text == decoded

            root.display_codes(freq)
            lav = sum(freq[c] * len(root.codebook[c]) for c in freq)
            lav = lav / sum(freq.values())
            print(f"Average code length: {lav}, Entropy:", entropy(freq))
