import pickle
import random
import os
import fileinput
import subprocess
import re
import sys
import base64

attempts = 0


def vec2string(vec):
    res = ""
    for e in vec:
        res += str(poly2int[e]) + ","
    return res[:-1]


def mat2strings(mat):
    res = []
    for row in mat:
        res.append(vec2string(row))
    return base64.b64encode("\n".join(res).encode()).decode()

def strings2vec(str):
    res = []
    vals = str.split(",")
    for val in vals:
        res.append(K(Integer(int(val)).digits(2)))
    return vector(K, res)


def strings2mat(strs, nrows, ncols):
    strs = base64.b64decode(strs.encode()).decode().split("\n")
    mat = []
    for row in strs:
        mat.append(strings2vec(row))
        assert ncols == len(mat[-1])
    assert nrows == len(mat)
    return matrix(K, mat)



assert len(sys.argv) >= 6, "q n m o2 PK0 PK1 PK2 ..."

# read PK
q, n, m, o2 = map(int, sys.argv[1:5])

K = GF(q)
poly2int = {K(Integer(i).digits(2)): i for i in range(q)}

PK = [strings2mat(s, n, n) for s in sys.argv[5:]]

if q != 16:
    print("Demo only for q = 16")
    exit()

# if (n-o2) % 2 != 0:
#    print("Demo only for parameters with n-o2 even")
#    exit()

if (n - o2) - 2 * (m - 2) > 2:
    print(
        "Demo only for when (n-o2)-2*(m-02) is not too large, otherwise Kipnis-Shamir is too slow."
    )
    exit()

# compile Wiedemann XL
print("Compiling Wiedemann XL ...")

with open("./xl-20160426/Makefile") as f:
    lines = f.readlines()

lines[0] = "Q = " + str(q) + "\n"
lines[1] = "M = " + str(m - 1) + "\n"
lines[2] = "N = " + str(n - m - 2) + "\n"

print("", lines[0], lines[1], lines[2])

with open("./xl-20160426/Makefile", "w") as f:
    f.writelines(lines)

subprocess.run("cd ./xl-20160426; make", shell=True)
subprocess.run("cp ./xl-20160426/xl .", shell=True)

load("./Rainbow.sage")
load("./SimpleAttack.sage")

while True:
    # make a guess and compose an MQ system
    print("Make a guess and compose MQ system")
    guess, system_filename, D_x_ker = Attack(PK)

    print("guess = ", guess)

    # Run Wiedemann XL
    print("Run Wiedemann XL")
    subprocess.run(
        "./xl --challenge " + system_filename + " --all | tee WXL_output.txt",
        shell=True,
    )

    matches = [line for line in open("WXL_output.txt") if re.search(r"is sol", line)]

    if len(matches) >= 1:
        break


solution_string = matches[0]

print("Solution found after %d attempts" % attempts)
print(solution_string)

print("Parsed solution:")

y = vector(
    [0]
    + [str_to_elt(str.lower(s[1:])) for s in solution_string.split()[0 : n - m - 2]]
    + [1]
)
print(y)

print("y")
print(D_x_ker.transpose() * y)
y = D_x_ker.transpose() * y

Eval_guess = Eval(PK, guess)
Eval_y = Eval(PK, y)

print(Eval_guess)
print(Eval_y)

alpha = Eval_y[0] / Eval_guess[0]
print("alpha:", alpha)
oil_vec = sqrt(alpha) * guess + y

print("vector in O2 is x = sqrt(alpha)*guess + y:")
print(oil_vec)
print("Sanity check: Pk(x) should be zero:")
print("Pk(x) = ", Eval(PK, oil_vec))


print("Finishing the attack:")

# Finding W space
basis_Fn = (K ** n).basis()
W_spanning_set = [Differential(PK, e, oil_vec) for e in basis_Fn]
W = matrix((K ** m).span(W_spanning_set).basis()).transpose()

# Finding O2 space
W_perp = matrix(W.kernel().basis())
P1 = [
    sum([W_perp[j, i] * PK[i] for i in range(m)]) for j in range(m - o2)
]  # P1 is inner layer of Rainbow

O2 = K ** n
for P in P1:
    O2 = O2.intersection((P + P.transpose()).kernel())

O2 = matrix(O2.basis()).transpose()

# Finding O1 space

# extend basis of O2 to a basis of K^n
Basis = O2
for e in basis_Fn:
    if e not in Basis.column_space():
        Basis = Basis.augment(matrix(K, n, 1, list(e)))

UOV_pk = [Basis.transpose() * P * Basis for P in P1]
for p in UOV_pk:
    Make_UD(p)
UOV_pk = [p[-n + o2 :, -n + o2 :] for p in UOV_pk]

print("Kipnis Shamir attack to find vectors in O1:")
# Kipnis-Shamir attack look for eigenvalues of MM untill we have a basis for Oil space
# We check if an eigenvalue is in the oil space by checking Eval(UOV_pk).is_zero()

OV_O_basis = []
while len(OV_O_basis) < m - o2:
    M = sum([K.random_element() * (P + P.transpose()) for P in UOV_pk])
    M0 = UOV_pk[0] + UOV_pk[0].transpose()

    to_delete = None
    if (n - o2) % 2 == 1:
        to_delete = random.randrange(n - o2)
        M = M.delete_rows([to_delete])
        M0 = M0.delete_rows([to_delete])
        M = M.delete_columns([to_delete])
        M0 = M0.delete_columns([to_delete])

    if not M.is_invertible():
        continue
    MM = M0 * M.inverse()

    cp = MM.characteristic_polynomial()
    for f, a in factor(cp):
        fMM_ker = f(MM).kernel()
        if fMM_ker.rank() == 2:
            b1, b2 = fMM_ker.basis()
            for v in K:
                o = b1 + v * b2
                if (n - o2) % 2 == 1:
                    o_list = list(o)
                    o_list.insert(to_delete, K(0))
                    o = vector(o_list)

                if Eval(UOV_pk, o).is_zero():
                    if not o in span(OV_O_basis, K):
                        OV_O_basis.append(o)
                        print("O1 vectors found: %d" % len(OV_O_basis))
                        break
                    else:
                        print("Vector is not new :(")
            if len(OV_O_basis) == m - o2:
                break
        if len(OV_O_basis) == m - o2:
            break


# extend basis to K^n
O1_basis = [Basis * vector([K(0)] * o2 + list(b)) for b in OV_O_basis] + O2.columns()
O1 = matrix(O1_basis).transpose()

print("Encoded recovered SK")
print(mat2strings(O1))
print(mat2strings(O2))
print(mat2strings(W))
print("------------------")
