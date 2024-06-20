import os
import sys
import marshal
import itertools
import argparse
from operator import itemgetter
from functools import partial
from collections import Counter
import heapq
import time

try:
    import cPickle as pickle
except:
    import pickle 

termchar = 17 # you can assume the byte 17 does not appear in the input file



#wrap all huffman stuff into a class
class HuffmanNode:
    def __init__(self, char, freq):
        self.char = char
        self.freq = freq
        self.left = None
        self.right = None

    def __lt__(self, other):
        return self.freq < other.freq

#build a min heap w the huffman node 
def build_huffman_tree(freq_dict):
    heap = [HuffmanNode(char, freq) for char, freq in freq_dict.items()]
    heapq.heapify(heap)

    while len(heap) > 1:
        left = heapq.heappop(heap)
        right = heapq.heappop(heap)
        merged = HuffmanNode(None, left.freq + right.freq)
        merged.left = left
        merged.right = right
        heapq.heappush(heap, merged) 

    return heap[0]

#build the codes 
def build_huffman_codes(node, current_code, codes):
    if node is None:
        return
    if node.char is not None:
        codes[node.char] = current_code
    build_huffman_codes(node.left, current_code + '0', codes)
    build_huffman_codes(node.right, current_code + '1', codes)

# This takes a sequence of bytes over which you can iterate, msg, 
# and returns a tuple (enc,\ ring) in which enc is the ASCII representation of the 
# Huffman-encoded message (e.g. "1001011") and ring is your ``decoder ring'' needed 
# to decompress that message.
def encode(msg):
    #print(f"encode msg {msg}")  
    freq_dict = Counter(msg)
    
    #call methods to build the huffman trees
    root = build_huffman_tree(freq_dict)
    codes = {} 
    build_huffman_codes(root, '', codes)

    encoded_msg = ""
    for char in msg:
        # check if char exists as a key in the codes dictionary
        if char in codes:
            encoded_msg += codes[char] 
    
    #set decoder ring to the dictionary 
    decoder_ring = codes
    #print(f"encoded_msg is {encoded_msg}\n decoder_ring is {decoder_ring}")

    return encoded_msg, decoder_ring


# This takes a string, cmsg, which must contain only 0s and 1s, and your 
# representation of the ``decoder ring'' ring, and returns a bytearray msg which 
# is the decompressed message. 

# This takes a sequence of bytes over which you can iterate, msg, and returns a tuple (compressed, ring) 
# in which compressed is a bytearray (containing the Huffman-coded message in binary, 
# and ring is again the ``decoder ring'' needed to decompress the message.
def compress(msg, useBWT):
    if useBWT:
        bwt_msg = bwt(msg)
        mtf_msg = mtf(bwt_msg)
        encoded_msg, decoder_ring = encode(mtf_msg)
    else:
        encoded_msg, decoder_ring = encode(msg)

    compressed = bytearray()
    byte_buffer = 0
    buffer_length = 0

    for bit in encoded_msg:
        byte_buffer = (byte_buffer << 1) | int(bit)
        buffer_length += 1

        if buffer_length == 8:
            # buffer is full, write it to the compressed output
            compressed.append(byte_buffer)
            byte_buffer = 0
            buffer_length = 0

    #if there are remaining bits in the buffer, write them to the compressed output
    if buffer_length > 0:
        byte_buffer <<= 8 - buffer_length
        compressed.append(byte_buffer)

    #print(compressed)
    return compressed, decoder_ring


def decode(cmsg, decoderRing):
    # Creates an array with the appropriate type so that the message can be decoded.
    decoded_bytes = bytearray()
    current_code = ""

    for bit in cmsg:
        current_code += bit
        # check if the current_code exists in the decoder_ring
        if current_code in decoderRing.values():
            #find the corresponding character for the code
            decoded_char = [char for char, code in decoderRing.items() if code == current_code][0]
            decoded_bytes.append(decoded_char)  #append the decoded character as a byte to the bytearray
            current_code = ""  #reset the current_code for the next code

    return decoded_bytes

def decompress(msg, decoderRing, useBWT):
    decoded_bits = ""
    decoded_msg = bytearray()
    

    #unshifting bits and appending them to decoded_bits
    for byte in msg:
        for i in range(8):
            decoded_bits = decoded_bits + str((byte >> 7) & 1)
            byte = byte << 1

    current_code = ""
    #iterate through the binary bits and try to find matching codes in the decoderRing
    for bit in decoded_bits:
        current_code += bit
        if current_code in decoderRing:
            decoded_msg.append(decoderRing[current_code])
            current_code = ""

    if useBWT:
        decoded_msg = mtf(decoded_msg)
        decoded_msg = ibwt(decoded_msg)

    return decoded_msg

# memory efficient iBWT
def ibwt(msg): #tried my best here i think this is causing the decompress to 0 not rlly sure
    if not msg:
        return bytearray() #basecase return null 

    #create a list of BWT characters
    bwt_list = list(msg)

    #sort the list to find the original order
    sorted_bwt = sorted(bwt_list)

    # Calculate the count of each character in the BWT
    char_count = {}
    for char in msg:
        if char not in char_count:
            char_count[char] = 0
        char_count[char] += 1

    #make a dictionary to map characters to their positions in the sorted BWT
    char_to_position = {}
    position = 0
    for char, count in char_count.items():
        for i in range(count):
            char_to_position[char] = (sorted_bwt.index(char, position), i)
            position = char_to_position[char][0] + 1

    #start with the first character in the BWT
    current_char = bwt_list[0]
    result = bytearray()

    #repeat until reconstructed the entire message
    for _ in range(len(msg)):
        result.append(current_char)
        current_char = bwt_list[char_to_position[current_char][0]]
    #print(result)
    
    return result


# Burrows-Wheeler Transform fncs
def radix_sort(values, key, step=0):
    sortedvals = []
    radix_stack = []
    radix_stack.append((values, key, step))

    while len(radix_stack) > 0:
        values, key, step = radix_stack.pop()
        if len(values) < 2:
            for value in values:
                sortedvals.append(value)
            continue

        bins = {}
        for value in values:
            bins.setdefault(key(value, step), []).append(value)

        for k in sorted(bins.keys()):
            radix_stack.append((bins[k], key, step + 1))
    return sortedvals
            
# memory efficient BWT
def bwt(msg):
    def bw_key(text, value, step):
        return text[(value + step) % len(text)]

    msg = msg + termchar.to_bytes(1, byteorder='big')

    bwtM = bytearray()

    rs = radix_sort(range(len(msg)), partial(bw_key, msg))
    for i in rs:
        bwtM.append(msg[i - 1])

    #print(bwtM[::-1]) 
    return bwtM[::-1]

# move-to-front encoding fncs
def mtf(msg):
    # Initialise the list of characters (i.e. the dictionary)
    dictionary = bytearray(range(256))
    
    # Transformation
    compressed_text = bytearray()
    rank = 0

    # read in each character
    for c in msg:
        rank = dictionary.index(c) # find the rank of the character in the dictionary
        compressed_text.append(rank) # update the encoded text
        
        # update the dictionary
        dictionary.pop(rank)
        dictionary.insert(0, c)

    #dictionary.sort() # sort dictionary
    return compressed_text # Return the encoded text as well as the dictionary

# inverse move-to-front
def imtf(compressed_msg):
    compressed_text = compressed_msg
    dictionary = bytearray(range(256))

    decompressed_img = bytearray()

    # read in each character of the encoded text
    for i in compressed_text:
        # read the rank of the character from dictionary
        decompressed_img.append(dictionary[i])
        
        # update dictionary
        e = dictionary[i]
        dictionary.pop(i)
        dictionary.insert(0, e)
        
    return decompressed_img


if __name__=='__main__':

    # argparse is an excellent library for parsing arguments to a python program
    parser = argparse.ArgumentParser(description='<Insert a cool name for your compression algorithm> compresses '
                                                 'binary and plain text files using the Burrows-Wheeler transform, '
                                                 'move-to-front coding, and Huffman coding.')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-c', action='store_true', help='Compresses a stream of bytes (e.g. file) into a bytes.')
    group.add_argument('-d', action='store_true', help='Decompresses a compressed file back into the original input')
    group.add_argument('-v', action='store_true', help='Encodes a stream of bytes (e.g. file) into a binary string'
                                                       ' using Huffman encoding.')
    group.add_argument('-w', action='store_true', help='Decodes a Huffman encoded binary string into bytes.')
    parser.add_argument('-i', '--input', help='Input file path', required=True)
    parser.add_argument('-o', '--output', help='Output file path', required=True)
    parser.add_argument('-b', '--binary', help='Use this option if the file is binary and therefore '
                                               'do not want to use the BWT.', action='store_true')

    args = parser.parse_args()

    compressing = args.c
    decompressing = args.d
    encoding = args.v
    decoding = args.w


    infile = args.input
    outfile = args.output
    useBWT = not args.binary

    assert os.path.exists(infile)

    if compressing or encoding:
        fp = open(infile, 'rb')
        sinput = fp.read()
        fp.close()
        if compressing:
            msg, tree = compress(sinput,useBWT)
            fcompressed = open(outfile, 'wb')
            marshal.dump((pickle.dumps(tree), msg), fcompressed)
            fcompressed.close()
        else:
            msg, tree = encode(sinput)
            print(msg)
            fcompressed = open(outfile, 'wb')
            marshal.dump((pickle.dumps(tree), msg), fcompressed)
            fcompressed.close()
    else:
        fp = open(infile, 'rb')
        pck, msg = marshal.load(fp)
        tree = pickle.loads(pck)
        fp.close()
        if decompressing:
            sinput = decompress(msg, tree, useBWT)
        else:
            sinput = decode(msg, tree)
            print(sinput)

        fp = open(outfile, 'wb')
        fp.write(sinput)
        fp.close()
