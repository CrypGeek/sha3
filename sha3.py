import numpy as np

#任意長のメッセージを固定長のハッシュ(十六進数)
#utf-8
#python3じゃないと動かないです

#ハッシュ値の文字数を設定
HashLength = 32
#入力メッセージ
Message    = "こんにちは、SHA-3です。"

SHA3_WORD = (1600/25)/8 #ワード長w SHA-3では64(bit)に設定されている
SHA3_BLOCKSIZE = 25*SHA3_WORD #ブロック長b 圧縮関数に入力するビットサイズ　SHAでは1600(bit)
SHA3_ROUND = int(12 + 2 * np.log2(SHA3_WORD * 8)) #ラウンドnr 攪拌関数内で何回攪拌するか
c = HashLength * 32 #capacity 衝突困難性の指標になる(ハッシュ長が長いと衝突しにくいからかと)
r = 1600 - c #メッセージを分割するbit数
MM1M5 = [4,0,1,2,3] #x-1 mod 5
MP1M5 = [1,2,3,4,0] #x+1 mod 5
MP2M5 = [2,3,4,0,1] #x+2 mod 5
M2XP3YM5 = [[0,2,4,1,3],[3,0,2,4,1],[1,3,0,2,4],[4,1,3,0,2],[2,4,1,3,0]] #2x + 3y mod 5
R  = [[0,  1, 62, 28, 27 ], [36, 44,  6, 55, 20 ], [3, 10, 43, 25, 39 ], [41, 45, 15, 21, 8 ], [18,  2, 61, 56, 14]] #定数
RC = [
	0x0000000000000001, 0x0000000000008082, 0x800000000000808A, 0x8000000080008000, 0x000000000000808B, 0x0000000080000001,
	0x8000000080008081, 0x8000000000008009, 0x000000000000008A, 0x0000000000000088, 0x0000000080008009, 0x000000008000000A,
	0x000000008000808B, 0x800000000000008B, 0x8000000000008089, 0x8000000000008003, 0x8000000000008002, 0x8000000000000080,
	0x000000000000800A, 0x800000008000000A, 0x8000000080008081, 0x8000000000008080, 0x0000000080000001, 0x8000000080008008
] #round constants

def ROT(val, r_bits):
    if(val == 0.0):
        return 0
    else:
        max_bits = int(np.log2(val)) + 1
        val = (val << r_bits%max_bits) & (2**max_bits-1) | ((val & (2**max_bits-1)) >> (max_bits-(r_bits%max_bits)))
        return val

def Round(MsgInt):
    A = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
    for i in range(25): #5×5の内部状態に分割 64bit毎
        A[i] = (MsgInt >> (1600 - 64*(i + 1))) & (2**64 -1)
    for rnd in range(SHA3_ROUND):
        B = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
        C = [0,0,0,0,0]
        D = [0,0,0,0,0]
        for x in range(5): #θ step
            C[x] = A[x + 0*5] ^ A[x + 1*5] ^ A[x + 2*5] ^ A[x + 3*5] ^ A[x + 4*5]
        for x in range(5):
            D[x] = C[MM1M5[x]] ^ ROT(C[MP1M5[x]],1)
        for x in range(5):
            for y in range(5):
                A[x + 5*y] = A[x + 5*y] ^ D[x]
        for x in range(5): #ρ and π step
            for y in range(5):
                B[y + 5*M2XP3YM5[x][y]] = ROT(A[x + 5*y],R[x][y])
        for x in range(5):#χ step
            for y in range(5):
                A[x + 5*y] = B[x + 5*y] ^ (~B[x + 1 + 5*y] & B[x + 2 + 5*y])
        A[0] ^= RC[rnd] #ι step
    return A

def split(S): #Sをa,bに分ける
    comS = 0
    for i in range(len(S)):
        comS = (comS << i*64) + S[i]
    an = comS >> r
    bn = comS & (2**r - 1)
    return an,bn

if(len(Message) > 64):
	print("メッセージは" + Message[:64] + "...")
else:
	print("メッセージは" + Message)

msgInt = int("0x"+ Message.encode().hex(),0) #文字列を数字に変換
pad = r - ((len(bin(msgInt)) - 2) % r) #パディングするbit数d
padMsgInt = msgInt << pad #padding
msgBitLen = len(bin(padMsgInt)) - 2 #bin(x) = "0b..."なので"0b"分引いてます
blockCount = int(msgBitLen / r) #攪拌関数にかける回数=メッセージの分割数

Block = []
for blocknum in range(blockCount): #メッセージをrごとに分割
    blockInt = (padMsgInt >> (msgBitLen - r*(blocknum + 1))) & (2**r -1)
    Block.append(blockInt)

m1XOR0 = Block[0] ^ 0
m1PIV = m1XOR0 << c #m1||IV 本当はIV||m1
S1 = Round(m1PIV)
a1,b1 = split(S1)
a = a1;b = b1

if(blockCount > 1):
    for block in Block[1:]:
        miXORb = block ^ b
        aPb = (a << r) + miXORb
        Si = Round(aPb)
        a,b = split(Si)

print("ハッシュ値は" + hex(a)[2:HashLength + 2])
