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

__all__ = ['PolyClass', 'VectorSet', 'Polynomial']


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

start_index = 1  # for display purposes only

# todo make all this stuff just tuples, maybe zip together entries/signs
# todo check the combination stuff (vector bases) and enumeration stuff


class Permutation:
    """ Object representing permutations. Permutations are automatically
    standardized to 0 through n-1 on creation, which is used extensively
    later on. """

    def __init__(self, entries, signs=None, b=None):  # TODO is b the right symbol

        if isinstance(entries, Permutation):  # clone constructor
            self.entries = entries.entries
            self.signs = entries.signs
            self.b = entries.b
            return

        assert len(set(entries)) == len(entries), 'make sure elements are distinct!'
        assert len(signs) == len(entries), 'length of signs and entries should match!'

        if isinstance(entries, tuple) or isinstance(entries, list):
            sorted_entries = sorted(entries)
            self.entries = tuple(map(lambda e: sorted_entries.index(e), entries))
        else:
            # TODO raise error
            pass

        if signs:
            self.signs = tuple(signs)
        else:
            # TODO raise error
            pass

        if b:
            self.b = b
        else:
            # TODO raise error
            pass

    def __repr__(self):
        return ''.join([str(a+start_index) + '\''*(b+1) + '\t' for a, b in zip(self.entries, self.signs)])

    def __eq__(self, other):
        if isinstance(other, Permutation):
            return self.entries == other.entries and self.signs == other.signs
        else:
            return False

    def __hash__(self):
        return hash(self.entries) + 2*hash(self.signs)

    def __len__(self):
        return len(self.entries)

    def delete(self, idx) -> 'Permutation':
        return Permutation(self.entries[:idx] + self.entries[idx+1:], self.signs[:idx] + self.signs[idx+1:], self.b)

    def plot(self):  # TODO update
        """ Draws a plot of the given Permutations. """
        n = len(self)
        array = [[' ' for i in range(n)] for j in range(n)]
        for i in range(n):
            array[self[i]][i] = '*'
        array.reverse()
        s = '\n'.join( (''.join(l) for l in array))
        # return s
        print(s)

    def fills(self):  # TODO check this
        pp = PegPerm(self, [0]*len(self))
        pp, vectorset = pp.clean()  # TODO this should call new method "reduce" that does what their "compact" did
        # the * unpacks the set, so we zip together vectors
        # the absolute value is magic: if there is a dot by itself, the compact
        # function puts a 0 in the vector, then the -1 and abs turns it to a 1
        vector = [abs(max(t)-1) for t in zip(*vectorset.basis)]
        # pp expanded by the vector should match the original permutation
        assert pp.inflate(vector) == self
        return pp, vector


class PegPerm:
    """ A permutation with 'pegged' entries, representing classes of permutations.
    Pegs can be +, ., or -. Object represented by a sign tuple and a permutation. """

    def __init__(self, perm, decorators=None):

        if isinstance(perm, PegPerm):  # clone constructor
            self.perm = Permutation(perm.perm)
            self.decorators = perm.decorators
            return

        if isinstance(perm, Permutation) and (isinstance(decorators, list) or isinstance(decorators, tuple)):
            self.perm = Permutation(perm)
        else:
            raise Exception('Bad inputs for PegPerm constructor')

        self.perm = perm
        self.decorators = tuple(decorators)  # TODO MAKE SURE THEY ARE +/- 1

    def __len__(self):
        return len(self.decorators)

    def __repr__(self):
        signsymbols = ['(.)', '(+)', '(-)']
        return ''.join([str(a+start_index) + '\''*(b+1) + signsymbols[c] + '\t'
                        for a, b, c in zip(self.perm.entries, self.perm.signs, self.decorators)])

    def __hash__(self):
        return hash(self.perm) + 2*hash(self.decorators)

    def __eq__(self, other):
        if isinstance(other, PegPerm):
            return self.decorators == other.decorators and self.perm == other.perm
        else:
            return False

    def plot(self):  # TODO FINISH
        n = self.__len__()
        array = [[' ' for i in range(n)] for j in range(n)]
        signlist = ['.', '+', '-']
        for i in range(n):
            array[self.perm[i]][i] = str(signlist[self.signs[i]])
        array.reverse()
        s = '\n'.join( (''.join(l) for l in array))
        # return s
        print(s)

    def minimum_fill(self):
        newentries = []
        newsigns = []

        for entry, sign, decorator in zip(self.perm.entries, self.perm.signs, self.decorators):
            newentries.append(3*entry)
            newsigns.append(sign)
            if decorator != 0:
                newentries.append(3*entry + decorator)
                newsigns.append(sign)

        return Permutation(newentries, newsigns)

    def inflate(self, vector) -> Permutation:
        assert len(vector) == len(self)

        new_entries = []
        new_signs = []

        multiplier = 2*max(vector)

        for entry, sign, decorator, number in zip(self.perm.entries, self.perm.signs, self.decorators, vector):
            for i in range(number):
                new_entries.append(multiplier*entry + decorator*i)
                new_signs.append(sign)
        return Permutation(new_entries, new_signs, self.perm.b)

    def delete(self, m):
        return PegPerm(self.perm.delete(m), self.decorators[:m]+self.decorators[m+1:])

    def dot(self, m):
        return PegPerm(self.perm, self.decorators[:m]+(0,)+self.decorators[m+1:])

    def dotall(self):
        p = self.perm
        n = len(self)
        S = set()

        for i in range(1, 2**n):

            # builds binary vector from i
            k = i
            v = [0 for x in range(n)]
            for m in range(n):
                if k >= 2**(n-m-1):
                    v[m]=1
                    k = k-2**(n-m-1)

            # dots permutation using the vector
            q = PegPerm(self)
            for j in range(n):
                if v[j] == 1:
                    q = q.dot(j)  # TODO this could certainly be faster if we construct all these perms 1 dot at a time

            # adds q to the set
            S.add(q)

        return S

    def complete(self) -> set['PegPerm']:  # TODO could this be faster?
        n = len(self)
        L = [set() for i in range(n + 1)]
        L[n] = self.dotall()
        L[n].add(self)
        while n >= 1:
            for q in L[n]:
                for i in range(len(q)):
                    if len(q) > 1:
                        L[n-1].add(q.delete(i))
            n -= 1
        return set().union(*L)

    def is_compact(self) -> bool:  # TODO write cleaner
        for i in range(len(self)-1):
            if self.perm.signs[i] == self.perm.signs[i+1]:
                a = (self.perm.entries[i]+2) * self.decorators[i]
                b = (self.perm.entries[i+1]+2) * self.decorators[i+1]

                if b-a == 1 or (b == 0 and self.perm.entries[i+1]-self.perm.entries[i] == self.decorators[i]) or (a == 0 and self.perm.entries[i+1]-self.perm.entries[i] == self.decorators[i+1]):
                    return False
        return True

    def clean(self) -> tuple['PegPerm', VectorSet]:  # intended to work on PegPerms that are known to be compact

        entries = []
        signs = []
        decorators = []
        vector = []  # essentially the "maximum" vector we could inflate by, where 0 represents infinity

        idx = 0
        while idx < len(self):
            if self.decorators[idx] == 0:
                sign = self.perm.signs[idx]

                runlength = 1
                while runlength+idx < len(self):
                    if self.perm.signs[idx+runlength] == sign:
                        if abs(self.perm.entries[idx+runlength] - self.perm.entries[idx+runlength-1]) == 1:
                            runlength += 1
                        else:
                            break
                    else:
                        break

                if runlength == 1:
                    entries.append(self.perm.entries[idx])
                    signs.append(self.perm.signs[idx])
                    decorators.append(self.decorators[idx])
                    vector.append(0)
                else:
                    entries.append(self.perm.entries[idx])
                    signs.append(sign)
                    decorators.append(self.perm.entries[idx+1]-self.perm.entries[idx])
                    vector.append(runlength)

                idx += runlength

            else:
                entries.append(self.perm.entries[idx])
                signs.append(self.perm.signs[idx])
                decorators.append(self.decorators[idx])
                vector.append(0)
                idx += 1

        # create avoidance basis
        basis = set()
        for idx, val in enumerate(vector):
            if val:
                new_vector = [0 for i in range(len(entries))]
                new_vector[idx] = val
                basis.add(tuple(new_vector))

        # minimum vector is just the minimum vector that fills our new perm
        # (idk if this is natural but that's what the paper describes)
        minimum = [1 + abs(i) for i in decorators]

        return PegPerm(Permutation(entries, signs, self.perm.b), decorators), VectorSet(minimum, basis)

    def split_blocks(self, *splits):  # helper
        splits = list(splits)
        splits.sort()
        splits.reverse()
        n = len(splits)

        new_decorators = list(self.decorators)
        mult = 3  # TODO this could prob be 1
        new_entries = [mult*n*x for x in self.perm.entries]
        new_signs = list(self.perm.signs)

        spacer = n
        for split in splits:
            sign, val = self.decorators[split], self.perm.entries[split]
            new_decorators.insert(split+1, sign)
            new_entries.insert(split+1, mult*n*val + spacer*sign)
            new_signs.insert(split+1, self.perm.signs[split])
            spacer -= 1
        return new_entries, new_signs, new_decorators

    def block_transpose(self):  # TODO later
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                for k in range(j, n):
                    if self.signs[k] == 0:
                        continue
                    signs, perm = self.split_blocks(i, j, k)
                    signs = signs[:i+1] + signs[j+2:k+3] + signs[i+1:j+2] + signs[k+3:]
                    perm = perm[:i+1] + perm[j+2:k+3] + perm[i+1:j+2] + perm[k+3:]
                    S.add(PegPerm(signs, perm))
        return S

    def cut_paste(self):  # TODO later
        ''' Returns the set of PegPerms which can result from one cut paste
                move from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                for k in range(j, n):
                    if self.signs[k] == 0:
                        continue
                    signs, perm = self.split_blocks(i, j, k)
                    # break up sign into pieces
                    s1 = signs[:i+1]; s2 = signs[j+2:k+3]
                    s3 = signs[i+1:j+2]; s4 = signs[k+3:]
                    # break up perm into pieces
                    p1 = perm[:i+1]; p2 = perm[j+2:k+3]
                    p3 = perm[i+1:j+2]; p4 = perm[k+3:]
                    S.add( PegPerm(s1 + s2 + s3 + s4, p1 + p2 + p3 + p4) )
                    # reverse piece two, add it to set
                    s2.reverse(); p2.reverse()
                    s2 = [-sign for sign in s2]
                    S.add( PegPerm(s1 + s2 + s3 + s4, p1 + p2 + p3 + p4) )
                    # undo the reversal, then do the same for piece 3
                    s2.reverse(); p2.reverse()
                    s2 = [-sign for sign in s2]
                    s3.reverse(); p3.reverse()
                    s3 = [-sign for sign in s3]
                    S.add( PegPerm(s1 + s2 + s3 + s4, p1 + p2 + p3 + p4) )
        return S

    def block_interchange(self):  # TODO later
        ''' Returns the set of PegPerms which can result from one block
                interchange from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                for k in range(j, n):
                    if self.signs[k] == 0:
                        continue
                    for m in range(k, n):
                        if self.signs[m] == 0:
                            continue
                        s, p = self.split_blocks(i, j, k, m)
                        s = s[:i+1] + s[k+3:m+4] + s[j+2:k+3] + s[i+1:j+2] + s[m+4:]
                        p = p[:i+1] + p[k+3:m+4] + p[j+2:k+3] + p[i+1:j+2] + p[m+4:]
                        S.add(PegPerm(s, p))
        return S

    def block_reversal(self):
        """ Returns the set of PegPerms which can result from one block
        reversal from the current one. """
        S = set([])
        n = len(self)

        for i in range(n):
            if self.decorators[i] == 0:
                continue
            for j in range(i, n):
                if self.decorators[j] == 0:
                    continue

                perm, signs, decorators = self.split_blocks(i, j)
                b = self.perm.b

                p_start, p_mid, p_end = perm[:i+1], perm[i+1:j+2], perm[j+2:]
                s_start, s_mid, s_end = signs[:i+1], [(s+1) % b for s in signs[i+1:j+2]], signs[j+2:]
                d_start, d_mid, d_end = decorators[:i+1], [-d for d in decorators[i+1:j+2]], decorators[j+2:]

                p_mid.reverse()
                s_mid.reverse()
                d_mid.reverse()

                S.add(PegPerm(Permutation(p_start + p_end, s_start + s_end, b), d_start + d_end))
        return S

    def prefix_transpose(self):  # TODO later
        ''' Returns the set of PegPerms which can result from one block
                reversal from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                signs, perm = self.split_blocks(i,j)
                s_start, s_mid, s_end = signs[:i+1], signs[i+1:j+2], signs[j+2:]
                p_start, p_mid, p_end = perm[:i+1], perm[i+1:j+2], perm[j+2:]
                signs = s_mid + s_start + s_end
                perm = p_mid + p_start + p_end
                S.add(PegPerm(signs, perm))
        return S

    def prefix_reversal(self):
        """ Returns the set of PegPerms which can result from one prefix
        reversal from the current one. """
        S = set()
        n = len(self)

        for i in range(n):
            if self.decorators[i] == 0:
                continue

            perm, signs, decorators = self.split_blocks(i)
            b = self.perm.b

            p_start, p_end = perm[:i+1], perm[i+1:]
            s_start, s_end = [(s+1) % b for s in signs[:i+1]], signs[i+1:]
            d_start, d_end = [-d for d in decorators[:i+1]], decorators[i+1:]

            p_start.reverse()
            s_start.reverse()
            d_start.reverse()

            S.add(PegPerm(Permutation(p_start+p_end, s_start+s_end, b), d_start+d_end))
        return S

    @staticmethod
    def sort_pr(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.prefix_reversal())
            S = S_
        return S

    @staticmethod
    def sort_pt(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.prefix_transpose())
            S = S_
        return S

    @staticmethod
    def sort_br(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.block_reversal())
            S = S_
        return S

    @staticmethod
    def sort_bt(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.block_transpose())
            S = S_
        return S

    @staticmethod
    def sort_bi(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.block_interchange())
            S = S_
        return S

    @staticmethod
    def sort_cp(k, b):
        S = {PegPerm(Permutation([0], [0], b), [1])}
        for i in range(k):
            S_ = set()
            for e in S:
                S_ = S_.union(e.cut_paste())
            S = S_
        return S

