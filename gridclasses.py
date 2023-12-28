#!/usr/bin/env python


import itertools as it
from math import factorial
from fractions import Fraction
# I'm sick of numpy, writing my own polynomial class with itertools
# from numpy import poly1d
import math
from operator import mul
from functools import reduce
import functools

from permutations import *


# Try to support both python 2 and 3 -
# some names are slightly different

import sys
if sys.version_info >= (3,0):
    izip_longest = it.zip_longest
else:
    izip_longest = it.izip_longest



using_sympy, using_sage = False, False


# check if we're running inside a sage session
using_sage = True
try: from sage.calculus.var import var
except: using_sage = False
finally: using_sympy = False

# if no sympy, check if we're in sage
if not using_sage:
    using_sympy = True
    try: from sympy import var, series, Matrix; using_sage = False
    except:
        using_sympy = False
        print('sympy not found, and not running in sage')
        print('generating functions require either sympy to be installed or this')
        print('module to be imported into a sage session')
        print('everything else will work fine\n')

# force it to print strings
# using_sage, using_sympy = False, False


# If not using sympy or sage, this variable decides whether or not to output
# long 'function strings', which can then be input into a CAS
# be warned, these can be VERY long
function_strings = True

# if true, then give exact bounds for when polynomial takes hold
# slows things down if you don't need it
exact_poly_bounds = False

# Changes the variable with which polynomials are displayed
poly_variable = 'n'


# if true, then give exact bounds for when polynomial takes hold
# slows things down if you don't need it
exact_poly_bounds = False


##########################################################
# Grid Class
##########################################################

class PolyClass(object):
    """ Represents a grid class as a dictionary mapping compact, clean PegPerms to VectorSets. """
    def __init__(self, arg1):

        # initialize with a cross-section dictionary
        if isinstance(arg1, dict):
            self.top_level_PegPerms = 'unknown'
            self.cross_sections = arg1

        # initialize with a single PegPerm
        # TODO: This is where performance is getting killed
        if isinstance(arg1, PegPerm):
            self.cross_sections = {}
            self.top_level = arg1
            PPSet = arg1.complete()
            for peg_perm in PPSet:
                if peg_perm.is_compact():
                    compact_perm, vector_set = peg_perm.clean()
                    if compact_perm in self.cross_sections:
                        self.cross_sections[compact_perm].union(vector_set)
                    else:
                        self.cross_sections[compact_perm] = vector_set

        # initialize with an iterable of PegPerms
        if isinstance(arg1, set) or isinstance(arg1, list):
            self.cross_sections = {}
            # for reference, this saves the original pegperms
            self.top_level = arg1
            for p in arg1:
                for peg_perm in p.complete():
                    if peg_perm.is_compact():
                        compact_perm, vector_set = peg_perm.clean()
                        if compact_perm in self.cross_sections:
                            self.cross_sections[compact_perm].union(vector_set)
                        else:
                            self.cross_sections[compact_perm] = vector_set

    def is_member(self, perm):
        pegperm, vector = perm.fills()
        return self.cross_sections[pegperm].is_member(vector)

    # these build grid classes from various block sorting methods
    @staticmethod
    def block_transpose(n, b):
        return PolyClass(PegPerm.sort_bt(n, b))

    @staticmethod
    def block_reversal(n, b):
        return PolyClass(PegPerm.sort_br(n, b))

    @staticmethod
    def prefix_reversal(n, b):
        return PolyClass(PegPerm.sort_pr(n, b))

    @staticmethod
    def prefix_transpose(n, b):
        return PolyClass(PegPerm.sort_pt(n, b))

    @staticmethod
    def block_interchange(n, b):
        return PolyClass(PegPerm.sort_bi(n, b))

    @staticmethod
    def cut_paste(n, b):
        return PolyClass(PegPerm.sort_cp(n, b))

    def polynomial(self):
        poly = 0
        start_term = 0
        for peg_perm in self.cross_sections:
            newpoly, start = self.cross_sections[peg_perm].polynomial()
            poly = poly + newpoly
            start_term = max(start_term, start)
        # if we have sympy or sage, we can figure out exactly when the polynomial
        # takes hold by comparing to the generating function terms
        if exact_poly_bounds and (using_sympy or using_sage):
            true_sequence = [term[1] for term in self.sequence(20)]
            poly_sequence = [poly(i) for i in range(1,21)]
            for i in range(1,20):
                if poly_sequence[i:] == true_sequence[i:]:
                    break
            start_term = i + 1
            print('polynomial takes hold at length %d' % start_term)
        else:
            print('polynomial works from at least length %d' % start_term)
        return poly

    def genfcn(self):
        if using_sympy or using_sage:
            fcn = 0
            for peg_perm, vectorset in self.cross_sections.items():
                # fcn = fcn + vectorset.generating_function().factor()
                fcn = fcn + vectorset.generating_function()
                # need to balance between memory usage and cpu usage
                # factoring often keeps memory free but uses cpu
            return fcn.factor()
        elif function_strings:
            biglist = []
            for peg_perm, vectorset in self.cross_sections.items():
                biglist += vectorset.generating_function()
            return biglist
        else:
            return 'sympy or sage required for symbolic generating fucntions'

    def sequence(self, n = 20):
        ''' If we can get a generating function, use that to get the exact sequence.
                otherwise, get a sequence from the polynomial which will be eventually
                correct. '''
        #TODO: when will the polynomial be correct?
        if using_sympy:
            gfcn = self.genfcn()
            # sympy magic, find series, extract coefficients
            coeffs = gfcn.series(x, 0, n + 1).as_poly().coeffs()
            coeffs.reverse()
            del coeffs[0]
            return list(enumerate(coeffs))[1:]
        if using_sage:
            gfcn = self.genfcn()
            # sage magic, find series, extract coefficients
            series = gfcn.series(x, n + 1)
            terms = series.coeffs(x)
            terms.sort(key = lambda term: term[1])
            coeffs = [term[0] for term in terms]
            return list(enumerate(coeffs))[1:]
        else:
            # compute the sequence from the polynomial, which will be innacurate for
            # the first ??? terms
            p = self.polynomial()
            return [int(p(i)) for i in range(1,n + 1)]

    def members(self, length):
        ''' Generates members of the polyclass by generating members of the
        vectorsets, then expanding the vectors into permutations. '''
        S = set()
        for pegperm, vectorset in self.cross_sections.items():
            for vector in vectorset.members(length):
                S.add(pegperm.expand(vector))
        return S


##########################################################
# VectorSet
##########################################################

class VectorSet():
    """ Represents the vector set as a minimum vector and a basis of avoidance vectors:
    {v | v >= minimum & v not >= b for all b in basis}. """

    def __init__(self, minimum, basis):
        assert isinstance(basis, set)
        for vector in basis:
            assert len(vector) == len(minimum)
            assert isinstance(vector, tuple)
            # assert all([ b >= m for b,m in zip(vector, minimum)])
        self.basis = basis
        self.minimum = tuple(minimum) # makes sure it's a tuple

    def __repr__(self):
        s = 'min = ' + self.minimum.__repr__() + '\n vector = '
        s += '\n'.join([b.__repr__() for b in self.basis])
        return s

    def __len__(self):
        ''' length of each vector (all have same length)'''
        return len(self.minimum)

    def basis_size(self):
        return len(self.basis)

    def is_member(self, vector):
        ''' returns true if vector is a member of the vector set '''
        if not all(v >= m for v,m in zip(vector, self.minimum)):
            return False
        for base in self.basis:
            if all(v >= b for v, b in zip(vector, base)):
                return False
        return True

    def members(self, length):
        ''' Generates members of the vectorset recursively. Should be able to speed
                it up... '''
        minsize = sum(self.minimum)
        if length < minsize:
            return set()
        if length == minsize:
            return set([self.minimum])
        else:
            S = set()
            for vector in self.members(length - 1):
                for i in (j for j in range(len(self.minimum)) if self.minimum[j] != 1):
                    newvec = list(vector)
                    newvec[i] += 1
                    if self.is_member(newvec):
                        S.add(tuple(newvec))
        return S

    def union(self, other):
        ''' Component-wise maximizes the basis elements. Modifies in place. '''
        newbasis = set([])
        for other_vector in other.basis:
            for my_vector in self.basis:
                # if either AVSet is empty, this is never called
                newbasis.add( tuple(map(max, zip(my_vector, other_vector))) )
        self.basis = newbasis
        self.clean_up()

    def clean_up(self):
        ''' Removes unnecessary basis elements '''
        newbasis = set()
        for vector in self.basis:
            # only want to keep a basis element if it avoids all others
            # make the new minimum just be zero, because we don't care if basis
            # vectors are above the original minimum
            others = VectorSet([0 for i in self.minimum], self.basis.difference([vector]))
            if others.is_member(vector):
                newbasis.add(vector)
        self.basis = newbasis

    def intersect(self, other_basis):
        ''' Adds in extra basis elements. '''
        if isinstance(other_basis, set):
            self.basis = self.basis.union(other_basis)
        elif isinstance(other_basis, tuple):
            self.basis.add(other_basis)
        else:
            print('other_basis should be a vector (tuple) or a set')

    def polynomial(self):
        ''' Builds a polynomial enumerating the class. '''
        k = len(self)
        dots = self.minimum.count(1)
        j = sum(self.minimum)
        poly = Polynomial.binom_poly(k - dots, j)
        basis_size = len(self.basis)
        mult = 1
        start_term = 0
        for level in range(1, basis_size + 1):
            mult *= -1
            for comb in it.combinations(self.basis, level):
                # comb is a list of vectors
                maxvector = [max(t) for t in zip(*comb)]
                maxvector = [max(t) for t in zip(maxvector, self.minimum)]
                poly = poly + mult * Polynomial.binom_poly(k - dots, sum(maxvector))
                start_term = max(start_term, sum(maxvector))
        return poly, start_term

    def generating_function(self):
        ''' starts with minimum, subtracts off basis '''
        if using_sage or using_sympy:
            k = len(self)
            dots = self.minimum.count(1)
            dimension = k - dots
            # sets up 'x' as an algebraic variable in both sympy and sage
            var('x')
            fcn = ( x ** sum(self.minimum) ) / ((1 - x) ** dimension)
            basis_size = len(self.basis)
            mult = 1
            for level in range(1, basis_size + 1):
                mult *= -1
                for comb in it.combinations(self.basis, level):
                    # comb is a list of vectors
                    maxvector = [max(t) for t in zip(*comb)]
                    maxvector = [max(t) for t in zip(maxvector, self.minimum)]
                    fcn = fcn + mult * (x ** sum(maxvector)) / ((1-x) ** dimension)
            return fcn
        if function_strings:
            # spits out LONG list of strings that can be input into sage
            # use L = self.generating_function(); print '\n +'.join(L)
            biglist = []
            k = len(self)
            dots = self.minimum.count(1)
            dimension = k - dots
            # sets up 'x' as an algebraic variable in both sympy and sage
            mult = 1
            biglist.append('(%d)*(x ** %d)/(1 - x) ** %d'
                                            % (mult, sum(self.minimum), dimension))
            basis_size = len(self.basis)
            for level in range(1, basis_size + 1):
                mult *= -1
                for comb in it.combinations(self.basis, level):
                    # comb is a list of vectors
                    maxvector = [max(t) for t in zip(*comb)]
                    maxvector = [max(t) for t in zip(maxvector, self.minimum)]
                    biglist.append('(%d)*(x ** %d)/(1 - x) ** %d'
                                                    % (mult, sum(maxvector), dimension))
            return biglist
        else:
            return '''either install sympy or run this script in a sage session for
                                generating functions'''


    def sequence(self, number_of_terms = 20):
        ''' Spits out the first few terms of the sequence enumerating the set. '''
        if using_sympy:
            gfcn = self.generating_function()
            # sympy magic, find series, extract coefficients
            coeffs = gfcn.series(x, 0, number_of_terms + 1).as_poly().coeffs()
            coeffs.reverse()
            return list(enumerate(coeffs))[1:]
        if using_sage:
            # TODO: this
            gfcn = self.genfcn()
            series = gfcn.series(x, n + 1)
            terms = series.coeffs(x)
            terms.sort(key = lambda term: term[1])
            coeffs = [term[0] for term in terms]
            return list(enumerate(coeffs))[1:]
        else:
            # first few terms may be inaccurate!
            poly = self.polynomial()
            return map(poly, range(1, number_of_terms))


##########################################################
# Polynomial
##########################################################

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
            s += str(coef) + poly_variable + '^' + str(exp + 1) + ' + '
        return s[:-3] # chops off the final ' + '

    def __call__(self, val):
        ''' allows polynomial to be called as a function '''
        result = 0
        for power, coef in enumerate(self):
            result += coef * ( val ** power )
        return int(result)

    def __add__(self, other):
        if isinstance(other, Polynomial):
            zipped = izip_longest(self, other, fillvalue=0)
            return Polynomial([a + b for a, b in zipped])
        else:
            newpoly = self[:]
            newpoly[0] += other
            return Polynomial(newpoly)

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
    def binom_poly(k,j):
        ''' returns the polynomial representing the coefficients of
                x^j / (1-x)^k. see poly.pdf for explanation'''
        if k < 1:
            return Polynomial([0])
        else: denominator = factorial(k-1)
        numerator = Polynomial([1])
        for i in range(k-1):
            numerator = numerator * Polynomial([- j + i + 1, 1])
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


# # this is for debugging, makes it easy to time functions
# import time
# class Timer(object):
#       def __enter__(self):
#           self.__start = time.time()
#       def __exit__(self, type, value, traceback):
#           self.__finish = time.time()
#           self.time = self.__finish - self.__start
#
# # if program is run alone, do...
# if __name__ == '__main__':
#       timer = Timer()
#       with timer:
#           A = PolyClass.block_interchange(2, timer = True)
#           print ''
#           L = A.genfcn()
#           print ''
#           f = open('genfcnstring.txt', 'w')
#           f.write(L[0] + '\\')
#           for line in L[1:]:
#               f.write('\n + ' + line + '\\')
#           f.close()
#
#       print timer.time



