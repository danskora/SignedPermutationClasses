#!/usr/bin/env python    #?????


import itertools as it
from math import factorial
from fractions import Fraction
import math
from operator import mul
from functools import reduce


import sys
if sys.version_info >= (3, 0):
    izip_longest = it.zip_longest
else:
    izip_longest = it.izip_longest


# check if we're running inside a sage session
try:
    from sage.calculus.var import var
except:
    using_sage = False
else:
    using_sage = True


# if no sage, check if sympy is installed
using_sympy = False
if not using_sage:
    try:
        from sympy import var, series, Matrix
    except:
        print('sympy not found, and not running in sage')
        print('generating functions require either sympy to be installed or this')
        print('module to be imported into a sage session')
        print('everything else will work fine\n')
    else:
        using_sympy = True


# If not using sympy or sage, this variable decides whether or not to output
# long 'function strings', which can then be input into a CAS
# be warned, these can be VERY long
FUNCTION_STRINGS = False

# Changes the variable with which polynomials are displayed
POLY_VARIABLE = 'n'


# GridClass object, represented as a set of compact peg perms
# ================================================================ #

class GridClass(object):
    ''' This class represents the grid class of a given set of SignedPerms'''

    def __init__(self, arg1):
        ''' We can define a grid class of a single permutation or of a set of
            permutations.  Completion and compacting are performed upon
            initialization'''

        # if given a SignedPerm
        if isinstance(arg1, SignedPerm):
            self.original_set = {arg1}
            self.new_set = set()
            for perm in arg1.complete():
                if perm.compact():
                    self.new_set.add(perm)

        # if given a set of SignedPerms
        elif isinstance(arg1, set):
            self.original_set = arg1
            self.new_set = set()
            for p in arg1:
                for perm in p.complete():
                    if perm.compact():
                        self.new_set.add(perm)

        # if given a list of Signed Perms
        elif isinstance(arg1, list):
            self.original_set = set(arg1)
            self.new_set = set()
            for p in arg1:
                for perm in p.complete():
                    if perm.compact():
                        self.new_set.add(perm)

    def is_member(self, perm):
        ''' Checks if a given SignedPerm is a member of the class. '''
        compact_perm, _ = perm.fills()
        return compact_perm in self.new_set

    # these build grid classes from various block sorting methods

    @staticmethod
    def block_reversal(n=1):
        return GridClass(SignedPerm.sort_br(n))

    @staticmethod
    def prefix_reversal(n=1):
        return GridClass(SignedPerm.sort_pr(n))

    @staticmethod
    def cut_paste(n=1):
        return GridClass(SignedPerm.sort_cp(n))

    def polynomial(self):
        ''' Builds a polynomial enumerating the class. '''
        poly = Polynomial([0])
        for perm in self.new_set:
            n = len(perm)
            poly = poly + Polynomial.binom_poly(n)
        return poly

    def genfcn(self):
        ''' Builds the generating function for the class '''
        if using_sympy or using_sage:
            fcn = 0
            # sets up 'x' as an algebraic variable in both sympy and sage
            var('x')
            for perm in self.new_set:
                n = len(perm)
                fcn = fcn + (x ** n) / ((1 - x) ** n)
            return fcn
        elif FUNCTION_STRINGS:
            # spits out LONG list of strings that can be input into sage
            # use L = self.generating_function(); print '\n +'.join(L)
            biglist = []
            for perm in self.new_set:
                n = len(perm)
                biglist.append('x ** %d / (1 - x) ** %d' % (n, n))
            return biglist
        else:
            return '''either install sympy or run this script in a sage session for
                                        generating functions'''

    def sequence(self, n=20):
        ''' Use the enumerating polynomial to get the first n terms in the sequence. '''
        p = self.polynomial()
        return [int(p(i)) for i in range(1, n + 1)]

    @staticmethod
    def find_vectors(length, weight):
        ''' Recursive method to generate all positive integer vectors of a specific length and weight.
            This is a helper method for members(). '''
        if length > weight:
            return []
        if length == weight:
            return [[1 for x in range(length)]]
        if length == 0:
            return []
        if length == 1:
            return [[weight]]

        vectors = []
        for i in range(1, weight - length + 2):
            for vector in GridClass.find_vectors(length - 1, weight - i):
                vectors.append([i]+vector)
        return vectors

    def members(self, length):
        """ Generates all members of the grid class with a certain length. """
        members = []
        for perm in self.new_set:
            for vector in GridClass.find_vectors(len(perm), length):
                members.append(perm.inflate(vector))
        return members

# End GridClass Class
#================================================================#


# SignedPerm object, represented as a pair of tuples
#================================================================#
class SignedPerm(object):
    """ A signed permutation object represented by a tuple of entries and a tuple of signs.
        The tuple of entries is automatically standardized to contain the numbers 1 through n. """
    def __init__(self, perm, signs, standardize=True):
        assert isinstance(perm, list) or isinstance(perm, tuple)
        assert isinstance(signs, list) or isinstance(signs, tuple)
        assert len(perm) == len(signs)
        self.signs = tuple(signs)
        for i in range(len(signs)):
            assert abs(signs[i]) == 1

        if standardize:
            self.perm = list(perm)
            assert len(set(perm)) == len(perm), 'make sure elements are distinct!'
            identity = sorted(self.perm)
            standardization = map(lambda e: identity.index(e) + 1, self.perm)
            self.perm = tuple(standardization)
        else:
            self.perm = tuple(perm)

    def __len__(self):
        """ Returns the length of the SignedPerm. """
        # sign vector and perm have the same length
        return len(self.signs)

    def __repr__(self):
        """ Changes to a string, for display. """
        if len(self) == 0:
            return 'e'
        signs = ['', '+', '-']
        string = ''.join([signs[y] + str(x) + ' ' for x,y in zip(self.perm, self.signs)])
        return string[:-1]

    def __hash__(self):
        """ Hashes a tuple containing perm tuple and signs tuple """
        return (self.signs, self.perm).__hash__()

    def __eq__(self, other):
        """ Determines whether two SignedPerms are equal. This method combined
        with the __hash__ method allows SignedPerms to be thrown into sets """
        return self.signs == other.signs and self.perm == other.perm

    def __getitem__(self, key):
        return self.perm[key] * self.signs[key]

    def plot(self):  #??????? Keep this or what
        """ Draws a plot of the given SignedPerm. """
        n = self.__len__()
        array = [[' ' for i in range(n)] for j in range(n)]
        signs = ['', '+', '-']
        for i in range(n):
            array[self.perm[i]][i] = str(signs[self.signs[i]])
        array.reverse()
        s = '\n'.join( (''.join(l) for l in array))
        print(s)

    def inflate(self, vector):
        """ Inflates each entry of the SignedPerm from a given integer vector. """

        assert len(vector) == len(self), ' vector length must match permutation length '
        if len(self) == 0:
            return SignedPerm([], [])
        new_perm = []
        new_signs = []
        assert min(vector) >= 0, ' vector cannot have negative entries '
        # we will capitalize on the autostandardization of SignedPerms
        spacer = 2*max(vector)
        for sign, entry, interval in zip(self.signs, self.perm, vector):
            for i in range(interval):
                new_perm.append(spacer * entry + sign * i)
                new_signs.append(sign)
        return SignedPerm(new_perm, new_signs)

    def delete(self, m):
        """ Returns the permutation which would result from removing index m.
                Does not modify the original permutation. """
        new_perm = list(self.perm)
        new_signs = list(self.signs)

        # let's standardize new_perm here since we can do it in O(n) rather than O(nlogn)
        entry = new_perm[m]
        del new_perm[m]
        del new_signs[m]
        for i in range(len(new_perm)):
            if new_perm[i] > entry:
                new_perm[i] -= 1

        return SignedPerm(new_perm, new_signs, standardize=False)

    def complete(self):
        """ Returns the set of SignedPerms contained by this one.  Works by deleting
                all 2^n possible subsets of the n entries """
        S = {SignedPerm(self.perm, self.signs, standardize=False)}
        for i in range(len(self)-1, -1, -1):
            S_prime = set()
            for perm in S:
                S_prime.add(perm.delete(i))
            S.update(S_prime)
        return S

    def compact(self):
        """ This method returns a boolean indicating whether the permutation is compact. """
        for i in range(len(self)-1):
            if self[i+1]-self[i] == 1:
                return False
        return True

    def fills(self):
        """ Returns the unique compact signed permutation which is filled by this one,
            as well as the vector we would inflate it by. """
        perm = SignedPerm(self.perm, self.signs, standardize=False)
        if len(perm) == 0:
            return perm, []
        vector = [1]
        i = 0
        while i < len(perm)-1:
            if perm[i+1]-perm[i] == 1:
                perm = perm.delete(i)
                vector[-1] += 1
            else:
                i += 1
                vector.append(1)
        assert perm.inflate(vector) == self
        return perm, vector

    def split_blocks(self, *splits):
        """ Helper method for performing block operations on SignedPerms.  Splits the SignedPerm into segments
            by cutting in the middle of the blocks indicated by *splits.  We do not convert these segments into
            SignedPerms, as this would lose information about each block's relative order. The splits variable
            should be comma separated splitting indices.  For example: p.split_blocks(0,0,3,5) """
        splits = list(splits)
        splits.sort()
        n = len(splits)

        # convert splits from a list of indices to a list of differences of indices
        splits.append(0)
        splits = [splits[i]-splits[i-1] for i in range(n)]

        spacer = n+1
        current_perm = [spacer*x for x in self.perm]
        current_signs = list(self.signs)
        perm_list = []
        signs_list = []

        for split in splits:
            perm_list.append(current_perm[:split+1])
            signs_list.append(current_signs[:split+1])

            current_perm = current_perm[split:]
            current_signs = current_signs[split:]
            current_perm[0] = current_perm[0]+current_signs[0]

        perm_list.append(current_perm)
        signs_list.append(current_signs)

        return perm_list, signs_list

    def prefix_reversal(self):
        ''' Returns the set of SignedPerms which can result from one block
                reversal from the current one. '''
        S = set()
        n = len(self)
        for i in range(n):
            perm_list, signs_list = self.split_blocks(i)
            S.add(SignedPerm(list(reversed(perm_list[0])) + perm_list[1],
                             list(reversed([-x for x in signs_list[0]])) + signs_list[1]))
        return S

    def block_reversal(self):
        """ Returns the set of SignedPerms which can result from one block
                reversal from the current one. """
        S = set()
        n = len(self)
        for i in range(n):
            for j in range(i, n):
                perm_list, signs_list = self.split_blocks(i, j)
                S.add(SignedPerm(perm_list[0] + list(reversed(perm_list[1])) + perm_list[2],
                                 signs_list[0] + list(reversed([-x for x in signs_list[1]])) + signs_list[2]))
        return S

    def cut_paste(self):
        """ Returns the set of SignedPerms which can result from one cut paste
                move from the current one. """
        S = set()
        n = len(self)
        for i in range(n):
            for j in range(i, n):
                for k in range(j, n):
                    perm_list, signs_list = self.split_blocks(i, j, k)
                    S.add(SignedPerm(perm_list[0] + perm_list[2] + perm_list[1] + perm_list[3],
                                     signs_list[0] + signs_list[2] + signs_list[1] + signs_list[3]))
                    S.add(SignedPerm(perm_list[0] + perm_list[2] + list(reversed(perm_list[1])) + perm_list[3],
                                     signs_list[0] + signs_list[2] + list(reversed([-x for x in signs_list[1]])) + signs_list[3]))
                    S.add(SignedPerm(perm_list[0] + list(reversed(perm_list[2])) + perm_list[1] + perm_list[3],
                                     signs_list[0] + list(reversed([-x for x in signs_list[2]])) + signs_list[1] + signs_list[3]))
        return S

    @staticmethod
    def sort_pr(n = 1):
        """ Returns set of SignedPerms representing n prefix reversals """
        L = [{SignedPerm([1], [1])}]
        for i in range(n):
            S = set()
            for p in L[i]:
                S = S.union(p.prefix_reversal())
            L.append(S)
        return L[n]

    @staticmethod
    def sort_br(n = 1):
        """ Returns set of SignedPerms representing n block reversals """
        L = [{SignedPerm([1], [1])}]
        for i in range(n):
            S = set()
            for p in L[i]:
                S = S.union(p.block_reversal())
            L.append(S)
        return L[n]

    @staticmethod
    def sort_cp(n = 1):
        """ Returns set of SignedPerms representing n cut-and-paste operations """
        L = [{SignedPerm([1], [1])}]
        for i in range(n):
            S = set()
            for p in L[i]:
                S = S.union(p.cut_paste())
            L.append(S)
        return L[n]

# End SignedPerm Class
#================================================================#

# This class was borrowed directly from Homberger and Vatter,
# who created it as a substitute for numpy
#================================================================#
class Polynomial(list):
    ''' Represents a polynomial as a list of coefficients.
            First coefficient is the constant term.  '''

    # stores this so it doesn't have to recompute multiple times
    base_change_matrix = None

    # uses regular list initialization, but trims off trailing zeroes
    def __init__(self, L):
        list.__init__(self, L)
        while self[-1] == 0 and len(self) > 1:
            del self[-1]
        if not self:
            self = [0]

    # displays as a polynomial
    def __repr__(self):
        s = str(self[0]) + ' + '
        for exp, coef in enumerate(self[1:]):
            s += str(coef) + POLY_VARIABLE + '^' + str(exp + 1) + ' + '
        return s[:-3] # chops off the final ' + '

    def __call__(self, val):
        ''' allows polynomial to be called as a function '''
        result = 0
        for power, coef in enumerate(self):
            result += coef * ( val ** power )
        return int(result)

    def __add__(self, other):
        if isinstance(other, Polynomial):
            zipped = izip_longest(self, other, fillvalue = 0)
            return Polynomial( [a + b for a,b in zipped] )
        else:
            newpoly = self[:]
            newpoly[0] += other
            return Polynomial( newpoly )

    def __mul__(self, other):
        if isinstance(other, Polynomial):
            cartesian_product = it.product(enumerate(self), enumerate(other))
            product = [0 for i in range(len(self) + len(other) - 1)]
            for first, second in cartesian_product: # magic
                product[first[0] + second[0]] += first[1] * second[1]
            return Polynomial(product)
        else:
            return Polynomial([other * coef for coef in self])

    @staticmethod
    def binom_poly(exp):
        ''' returns the polynomial representing the coefficients of
                (x)^exp / (1-x)^exp '''
        if exp < 1:
            return Polynomial([0])
        denominator = factorial(exp-1)
        numerator = Polynomial([1])
        for i in range(exp-1):
            numerator = numerator * Polynomial([- exp + i + 1, 1])
        return Fraction(1, denominator) * numerator

    @staticmethod
    def choose(k):
        ''' Returns the polynomial (n choose k) '''
        if k == 0:
            return Polynomial([1])
        else:
            numerator = reduce( mul, [Polynomial([-i, 1]) for i in range(k)])
            return Fraction(1, factorial(k)) * numerator

    @staticmethod
    def generate_base_change_matrix(size = 8):
        ''' Builds the matrix. '''
        if using_sympy:
            A = Matrix([ Polynomial.choose(k).pad_with_zeroes(size)
                                                            for k in range(size) ])
            A = A.transpose()
            Polynomial.base_change_matrix = A.inv()
        elif using_sage:
            print('eh, ill do this later')
        else:
            print('need sympy or sage')


    def binom_basis(self):
        ''' Changes basis from {1, n, n^2 _ to (n choose 0), (n choose 1), ... '''
        if not Polynomial.base_change_matrix:
            Polynomial.generate_base_change_matrix(len(self))
        size = len(Polynomial.base_change_matrix)
        if len(self)**2 > size:
            Polynomial.generate_base_change_matrix(len(self))
        if not using_sympy:
            return 'use sympy'
        v = Matrix(self.pad_with_zeroes(int(math.sqrt(size))))
        return Polynomial.base_change_matrix * v


    def pad_with_zeroes(self, length):
        return list(self) + [0] * (length - len(self))


    # These all make the object work with built in operations, like adding,
    # subtracting, raising to powers, mulitplication, and assignment.

    def __rmul__(self, other):
        return self * other

    def __neg__(self):
        return -1 * self

    def __pow__(self, power):
        result = 1
        for i in range(power):
            result = result * self
        return result

    def __sub__(self, other):
        return self + (-1 * other)

    def __radd__(self, other):
        return self + other

    def __rsub__(self, other):
        return self - other
# End Polynomial Class
#================================================================#