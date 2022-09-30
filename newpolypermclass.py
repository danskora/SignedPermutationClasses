#!/usr/bin/env python


import itertools as it
from math import factorial
from fractions import Fraction
import math
from operator import mul
from functools import reduce



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

# Changes whether displayed permutations start at 1 or 0.
# Only affects the display, internally permutations always start at 0.
start_index = 1

# Changes the variable with which polynomials are displayed
poly_variable = 'n'


#i f true, then give exact bounds for when polynomial takes hold
# slows things down if you don't need it
exact_poly_bounds = False

# PolyClass object, represented as a dict of compact PegPerms and VectorSets
#================================================================#
class PolyClass(object):
    ''' Represents a polynomial class as a union of cross sections'''

    # TODO: add lots of ways of initializing (this will be a big function)
    def __init__(self, arg1, timer = False):

        # if given a dictionary, assume it's the cross section dictionary
        if isinstance(arg1, dict):
            self.top_level_PegPerms = 'unknown'
            self.cross_sections = arg1

        # if given a PegPerm, complete it and compact it first
        # TODO: This is where performance is getting killed
        if isinstance(arg1, PegPerm):
            self.cross_sections = {}
            self.top_level = arg1
            PPSet = arg1.complete()
            for peg_perm in PPSet:
                compact_perm, vector_set = peg_perm.compact()
                if compact_perm in self.cross_sections:
                    self.cross_sections[compact_perm].union(vector_set)
                else:
                    self.cross_sections[compact_perm] = vector_set

        # if given a set or list or pegperms, complete and compact them all
        if isinstance(arg1, set) or isinstance(arg1, list):
            if timer: i = 1
            self.cross_sections = {}
            # for reference, this saves the original pegperms
            self.top_level = arg1
            for e in arg1:
                print(e)
            for p in arg1:
                if timer: print(i); i += 1

                for peg_perm in p.complete():
                    compact_perm, vector_set = peg_perm.compact()
                    if compact_perm in self.cross_sections:
                        self.cross_sections[compact_perm].union(vector_set)
                    else:
                        self.cross_sections[compact_perm] = vector_set
            one = 0
            two = 0
            three = 0
            four = 0
            five = 0
            six = 0
            seven = 0
            eight = 0


            for key in self.cross_sections:
                match len(key):
                    case 1:
                        one+=1
                    case 2:
                        two+=1
                    case 3:
                        three+=1
                    case 4:
                        four+=1
                    case 5:
                        five+=1
                    case 6:
                        six+=1
                    case 7:
                        seven+=1
                    case 8:
                        eight+=1
            print("8: " + str(eight))
            print("7: " + str(seven))
            print("6: " + str(six))
            print("5: " + str(five))
            print("4: " + str(four))
            print("3: " + str(three))
            print("2: " + str(two))
            print("1: " + str(one))



            # with_base = 0
            # without_base = 0
            # for key in self.cross_sections:
            #     if len(key) != 0:
            #         if len(key) == list(key.perm)[-1] + 1 and 1 == list(key.signs)[-1]:
            #             with_base += 1
            #         else:
            #             without_base += 1
            #     else:
            #         without_base += 1
            #
            # print("With Base: " + str(with_base))
            # print("Without Base: " + str(without_base))


    def is_member(self, perm):
        ''' Checks if a given permutation is a member of the class. '''
        pegperm, vector = perm.fills()
        return self.cross_sections[pegperm].is_member(vector)



    # these build polyclasses from various block sorting methods
    # setting timer = True shows the progress
    # TODO: maybe too much, but a progress bar would be cute
    @staticmethod
    def block_transpose(n = 1, timer = False):
        return PolyClass(PegPerm.sort_bt(n), timer)

    @staticmethod
    def block_reversal(n = 1, timer = False):
        return PolyClass(PegPerm.sort_br(n), timer)

    @staticmethod
    def prefix_reversal(n = 1, timer = False):
        return PolyClass(PegPerm.sort_pr(n), timer)

    @staticmethod
    def prefix_transpose(n = 1, timer = False):
        return PolyClass(PegPerm.sort_pt(n), timer)

    @staticmethod
    def block_interchange(n = 1, timer = False):
        return PolyClass(PegPerm.sort_bi(n), timer)

    @staticmethod
    def cut_paste(n = 1, timer = False):
        return PolyClass(PegPerm.sort_cp(n), timer)

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

# End PolyClass Class
#================================================================#


# Permutation object, an extension of built-in tuple
#================================================================#
class Permutation(tuple):
    ''' Object representing permutations, built on python tuples. Permutations are
            automatically standardized to 0 through n-1 on creation, which is used
            extensively later on.'''

    # use __new__ instead of __init__ because tuples are immutable, so we must
    # create a new one rather than modify an existing one.
    # note that tuples come with their own __hash__ function and __eq__ function
    def __new__(cls, p, n = None):
        ''' Initializes a permutation object, internal indexing starts at zero. '''
        if isinstance(p, Permutation):
            return tuple.__new__(cls, p)
        elif isinstance(p, tuple):
            entries = list(p)[:]
        elif isinstance(p, list):
            entries = p[:]
        # standardizes, starting at zero
        assert len(set(entries)) == len(entries), 'make sure elements are distinct!'
        entries.sort()
        standardization =  map(lambda e: entries.index(e), p)
        return tuple.__new__(cls, standardization)

    def __repr__(self):
        ''' Tells python how to display a permutation object. Uses start_index
                variable. Note that internal indexing always starts at zero. '''
        return ''.join([str(i + start_index) + ' ' for i in self])

    # indexing starts at zero, be careful!
    def delete(self, index):
        ''' Returns a new permutation with the 'index' entry deleted. '''
        L = list(self)
        entry = L[index]
        del L[index]
        # uses the auto-standardization of the constructor
        return Permutation(L)

    def plot(self):
        ''' Draws a plot of the given Permutations. '''
        n = self.__len__()
        array = [[' ' for i in range(n)] for j in range(n)]
        for i in range(n):
            array[self[i]][i] = '*'
        array.reverse()
        s = '\n'.join( (''.join(l) for l in array))
        # return s
        print(s)

    def fills(self):
        ''' Returns the PegPerm representing the class which is filled by this
                permutation. '''
        pp = PegPerm([0 for i in self], self)
        pp, vectorset = pp.compact()
        # the * unpacks the set, so we zip together vectors
        # the absolute value is magic: if there is a dot by itself, the compact
        # function puts a 0 in the vector, then the -1 and abs turns it to a 1
        vector = [abs(max(t)-1) for t in zip(*vectorset.basis)]
        # pp expanded by the vector should match the original permutation
        assert pp.expand(vector) == self
        return pp, vector


# End Permutation Class
#================================================================#




# PegPerm object, represented as a pair of tuples
#================================================================#
class PegPerm(object):

    ''' A permutation with 'pegged' entries, representing classes of permutations.
            Pegs can be +,., or -. Object represented by a sign tuple and a
            permutation '''
    def __init__(self, signs, perm = None):
        # if signs is actually already a PegPerm, just copy it
        if isinstance(signs, PegPerm):
            self.perm = Permutation(signs.perm)
            self.signs = signs.signs[:] # force copy by value, just in case
        elif perm == None:
            # tries to catch old code...
            # TODO: delete this
            raise Exception('update PegPerm initializer!')
        else:
            if isinstance(perm, list) or isinstance(perm, tuple):
                perm = Permutation(perm)
            else:
                raise Exception('PegPerm(sign, list)')
            assert isinstance(perm, Permutation)
            assert isinstance(signs, list) or isinstance(signs, tuple)
            assert len(signs) == len(perm)
            self.perm = perm
            self.signs = tuple(signs)

    def __len__(self):
        ''' Returns the length of the PegPerm. '''
        # sign vector and perm have the same length
        return len(self.signs)

    def __repr__(self):
        ''' Changes to a string, for display. '''
        signsymbols = [ '.', '+', '-' ]
        # list[-1] returns the last element of a list
        return ''.join( [signsymbols[x] + str(y + start_index)
                                        for x,y in zip(self.signs, self.perm)] )

    def __hash__(self):
        ''' Hashes a tuple containing signs tuple, perm tuple, and length. '''
        return (self.signs, self.perm).__hash__()


    def __eq__(self, other):
        ''' Determines whether two PegPerms are equal. This method combined
        with the __hash__ method allows PegPerms to be thrown into sets '''
        # broke this up, to hopefully make it faster
        return self.signs == other.signs and self.perm == other.perm
        # return [self.signs, self.perm] == [other.signs, other.perm]

    def plot(self):
        ''' Draws a plot of the given PegPerm. '''
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
        ''' Outputs the smallest permutation which is contained in this class but not
                contained in any smaller class. '''
        newperm = []
        for sign, entry in zip(self.signs, self.perm):
            newperm.add(3 * entry) # expand by 3 to avoid collisions
            if sign == 0:
                # do nothing if sign is a dot
                pass
            else:
                # expand entry into an increasing or decreasing bond
                newperm.add(3 * entry + sign)
        # the Permutation initializer standardizes automatically
        return Permutation(newperm)

    def expand(self, vector = None):
        ''' Expands each entry of the PegPerm from a given integer vector. '''
        # if no vector given, just return the smallest filling permutation
        if not vector:
            return self.minimum_fill()
        newperm = []
        # expands to avoid collisions
        multiplier = 2 * max(vector)
        assert len(vector) == len(self)
        for sign, entry, number in zip(self.signs, self.perm, vector):
            if sign == 0:
                assert number <= 1, 'trying to fit too much into a dot'
            for i in range(number):
                newperm.append(multiplier * entry + sign * i)
        return Permutation(newperm)

    def delete(self,m):
        ''' Deletes an entry of the permutation, and renumbers accordingly.
                Returns a new permutation.'''
        newperm = list(self.perm)
        newsign = list(self.signs)
        del newperm[m]
        del newsign[m]
        return PegPerm(newsign,newperm)

    def dot(self,m):
        ''' Changes any sign into a dot, returns new permutation. '''
        signs = list(self.signs)
        perm = self.perm
        signs[m] = 0
        return PegPerm(signs,perm)

    def dotall(self):
        ''' Changes signs to dots in all possible ways by building every 0,1 vector
        and dotting each 1 entry returns a SetOfTuples object. '''
        p = self.perm
        n = len(self)
        S = set([])
        # builds all possible 0,1 vectors, dots each '1' entry
        # uses binary decomposition of integers
        for i in range(1,2**n ):
            k = i
            v = [0 for i in range(n)]
            # builds binary vector from i
            for m in range(n):
                if k >= 2**(n-m-1):
                    v[m]=1
                    k = k-2**(n-m-1)
            # dots permutation using the vector
            q = PegPerm(self)
            for j in range(n):
                if v[j] ==1:
                    q = q.dot(j)
            S.add(q)
        return S

    def complete(self):
        ''' Uses deletions and (signs -> dots) to make a complete downset. Returns a
                PegPermSet '''
        # works by first dotting in all possible ways, then deleting entries from
        # each of the resulting PegPerms
        n = len(self)
        L = [set() for i in range(n + 1)]
        L[n].add(self)
        while n >= 1:
            for q in L[n]:
                for i in range(len(q)):
                    L[n-1].add(q.delete(i))          # NOTE Could be faster if we simply deleted all 2^n subsets of pegs rather than doing this thing recursively
            n -= 1
        return set().union(*L)

    def compact(self):
        ''' Returns an ElemSet. changes runs of dots to signs, absorbs dots into
        signs, and builds avoidance vector '''
        # This is where the work is done. It could be made shorter, but at a loss of
        # clarity. Works in stages: First it finds runs of dots and shrinks them
        # into signs ("cleaning" the peg permutation), and builds an avoidance vector.
        # Next it absorbs dots into adjacent signs (making the peg permutation compact),
        # if possible. Finally it splits up the avoidance vector.
        perm = self.perm
        # change to a list, to allow reassignment
        signs = list(self.signs)                 # NOTE remove cleanng step.  remove minimum vector bs
        vector = [0 for i in range(len(perm))]
        idx = 0
        # check for runs of dots, delete them and add to vector
        while idx < len(perm) - 1:
            val = perm[idx + 1] - perm[idx]
            # if a bond, checks if both signs are dots. if not, continues
            if abs(val) == 1:
                runlength = 1
                while (signs[idx + runlength - 1] == 0 and
                        signs[idx + runlength] == 0 and
                        perm[idx + runlength] - perm[idx + runlength - 1] == val):
                    runlength += 1
                    # stop if we hit the end of the permutation
                    if idx + runlength >= len(signs):
                        break
                # if we have a run, delete those entries and add to vector
                if runlength > 1:
                    for i in range(runlength - 1):
                        perm = perm.delete(idx + 1)
                        del signs[idx + 1]
                    del vector[idx + 1 : idx + runlength]
                    vector[idx] = runlength + 1
                    signs[idx] = val
            idx += 1
        # stage one finished
        # checks for off by one errors
        assert len(vector) == len(perm)
        assert len(perm) == len(signs)
        # now it merges stray dots into signs, if possible
        idx = 0
        while idx < len(perm) - 1:
            val = perm[idx + 1] - perm[idx]
            if val == 1:
                if ((signs[idx] == 1 and signs[idx + 1] != -1) or
                         (signs[idx] != -1 and signs[idx + 1] == 1)):
                    perm = perm.delete(idx + 1)
                    del signs[idx + 1]
                    signs[idx] = 1
                    if vector[idx + 1] * vector[idx] == 0:
                        del vector[idx + 1]
                        vector[idx] = 0
                    else:
                        newval = vector[idx] + vector[idx + 1]
                        del vector[idx + 1]
                        vector[idx] = newval
                else: idx += 1
            elif val == -1:
                if ((signs[idx] == -1 and signs[idx + 1] != 1) or
                         (signs[idx] != 1 and signs[idx + 1] == -1)):
                    perm = perm.delete(idx + 1)
                    del signs[idx + 1]
                    signs[idx] = -1
                    if vector[idx + 1] * vector[idx] == 0:
                        del vector[idx + 1]
                        vector[idx] = 0
                    else:
                        newval = vector[idx] + vector[idx + 1]
                        del vector[idx + 1]
                        vector[idx] = newval
                else: idx += 1
            else:
                idx += 1
        # at this point, we have a compact pegged permutation with a vector basis
        # need to split up vector into multiples if necessary
        basis = set()
        for idx, val in enumerate(vector):
            if val:
                basis_vector = [0 for i in range(len(signs))]
                basis_vector[idx] = val
                basis.add(tuple(basis_vector))
        # this makes a minimum vector
        minimum = [1 for i in signs]
        return PegPerm(signs, perm), VectorSet(minimum, basis)

    def split_blocks(self, *splits):
        ''' Splits some blocks into two pieces, used for block sorting.  The splits
                variable should be comma separated splitting points.
                For example: p.split_blocks(1,1,4,5) '''
        splits = list(splits)
        assert all( [self.signs[split] != 0 for split in splits] ), \
                                            "can't split dots!"
        # split bigger entries first to make indexing easier
        splits.sort(); splits.reverse()
        n = len(splits)
        newsigns = list(self.signs)
        # space out perm entries to avoid collisions
        # three is probably overkill, but no harm in making it big
        mult = 3
        newperm = [mult*n*x for x in self.perm]
        spacer = n
        for split in splits:
            sign, val = self.signs[split], self.perm[split]
            # don't allow splitting a dot
            newsigns.insert(split + 1, sign)
            newperm.insert(split + 1, mult*n*val + spacer*sign)
            spacer -= 1
        return newsigns, newperm

    def block_transpose(self):
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
                    signs, perm = self.split_blocks(i, j, k)
                    signs = signs[:i+1] + signs[j+2:k+3] + signs[i+1:j+2] + signs[k+3:]
                    perm = perm[:i+1] + perm[j+2:k+3] + perm[i+1:j+2] + perm[k+3:]
                    S.add(PegPerm(signs, perm))
        return S

    def cut_paste(self):
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

    def block_interchange(self):
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
        ''' Returns the set of PegPerms which can result from one block
                reversal from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            for j in range(i, n):
                if self.signs[j] == 0:
                    continue
                signs, perm = self.split_blocks(i,j)
                s_start, s_mid, s_end = signs[:i+1], signs[i+1:j+2], signs[j+2:]
                s_mid = [-sign for sign in s_mid[::-1]]
                signs = s_start + s_mid + s_end
                p_start, p_mid, p_end = perm[:i+1], perm[i+1:j+2], perm[j+2:]
                p_mid.reverse()
                perm = p_start + p_mid + p_end
                S.add(PegPerm(signs, perm))
        return S

    def prefix_transpose(self):
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
        ''' Returns the set of PegPerms which can result from one block
                reversal from the current one. '''
        S = set([])
        n = len(self)
        for i in range(n):
            if self.signs[i] == 0:
                continue
            signs, perm = self.split_blocks(i)
            s_start, s_end = signs[:i+1], signs[i+1:]
            p_start, p_end = perm[:i+1], perm[i+1:]
            p_start.reverse(), s_start.reverse()
            s_start = [-sign for sign in s_start]
            S.add(PegPerm(s_start + s_end, p_start + p_end))
        return S


    @staticmethod
    def sort_bt(n = 1):
        ''' Returns set of pegperms which result from n block transpositions from
                the identity. '''
        L = [set([PegPerm([1],[0])])]
        S = set([])
        for i in range(n):
            S = set([])
            for p in L[i]:
                S = S.union(p.block_transpose())
            L.append(S)
        return S

    @staticmethod
    def sort_pr(n = 1):
        ''' Returns set of pegperms which result from n prefix reversals from
                the identity. '''
        L = [set([PegPerm([1],[0])])]
        S = set([])
        for i in range(n):
            S = set([])
            for p in L[i]:
                S = S.union(p.prefix_reversal())
            L.append(S)
        return S

    @staticmethod
    def sort_pt(n = 1):
        ''' Returns set of pegperms which result from n prefix reversals from
                the identity. '''
        L = [set([PegPerm([1],[0])])]
        S = set([])
        for i in range(n):
            S = set([])
            for p in L[i]:
                S = S.union(p.prefix_transpose())
            L.append(S)
        return S

    @staticmethod
    def sort_br(n = 1):
        ''' Returns set of pegperms which result from n block reversals from
                the identity. '''
        L = [set([PegPerm([1],[0])])]
        S = set([])
        for i in range(n):
            S = set([])
            for p in L[i]:
                S = S.union(p.block_reversal())
            L.append(S)
        return S

    @staticmethod
    def sort_cp(n = 1):
        ''' Returns set of pegperms which result from n cut paste moves from
                the identity. '''
        L = [set([PegPerm([1],[0])])]
        S = set([])
        for i in range(n):
            S = set([])
            for p in L[i]:
                S = S.union(p.cut_paste())
            L.append(S)
        return S

    @staticmethod
    def sort_bi(n = 1):
        ''' Returns set of pegperms which result from n block interchanges from
                the identity. '''
        L = [set([PegPerm([1],[0])])]
        S = set([])
        for i in range(n):
            S = set([])
            for p in L[i]:
                S = S.union(p.block_interchange())
            L.append(S)
        return S
# End PegPerm Class
#================================================================#



# VectorSet object, represented as a minimum vector and a basis
#================================================================#
class VectorSet(object):
    ''' Represents the vector set
            {v | v >= minimum & v >\= basis for all vectors in basis} '''

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
        dots = 0
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
            dots = 0
            dimension = k - dots
            # sets up 'x' as an algebraic variable in both sympy and sage
            var('x')
            fcn = ( x ** sum(self.minimum) ) / ((1 - x) ** dimension)  #THIS IS WHERE I LEFT OFF
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
# End VectorSet Class
#================================================================#


# numpy likes to change integers to floats, seemingly at random.
# It drove me crazy so I wrote my own simple polynomial class.
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











