import copy
import gurobipy as gp


class Variable:
    created_variables = {}
    const = None

    def __init__(self, name, vartype):
        self.name = name
        self.type = vartype
        assert self.name not in Variable.created_variables.keys()
        Variable.created_variables[self.name] = self.type

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return self.name

    def __hash__(self):
        return hash(self.name)


Variable.const = Variable("ONE", 'c')


class Monomial:
    def __init__(self, variable, coefficient):
        self.variable = variable
        self.coefficient = coefficient

    def __add__(self, other):
        if self.alike(other):
            return Monomial(self.variable, self.coefficient + other.coefficient)
        assert False

    def scale(self, scalar):
        return Monomial(self.variable, self.coefficient * scalar)

    def alike(self, other):
        return self.variable == other.variable

    def __copy__(self):
        return Monomial(self.variable, self.coefficient)

    def __str__(self):
        if str(self.variable) != "":
            return f"{self.coefficient} {self.variable}"
        else:
            return str(self.coefficient)


class LinearCombination:
    def __init__(self, monomials):
        self.monomials = {}
        for m in monomials:
            assert isinstance(m, Monomial)
            self.monomials[m.variable] = m

    def __add__(self, other):
        new_vars = copy.copy(other)
        for var in self.monomials.keys():
            if new_vars.monomials.get(var, None) is not None:
                new_vars.monomials[var] += self.monomials[var]
            else:
                new_vars.monomials[var] = self.monomials[var]
        return new_vars

    def scale(self, scalar):
        new_vars = []
        for var in self.monomials.values():
            new_vars.append(var.scale(scalar))
        return LinearCombination(new_vars)

    def __mod__(self, other):
        result = []
        for m in self.monomials.values():
            result.append(Monomial(m.variable, m.coefficient % other))
        return LinearCombination(result)

    def __str__(self):
        s = ""
        for m in self.monomials.values():
            if m.coefficient != 0:
                s += str(m)
                s += " + "
        s = s.replace("+ -", "- ")
        return s[:-3]

    def __copy__(self):
        result = []
        for m in self.monomials.values():
            result.append(copy.copy(m))
        return LinearCombination(result)


class FormalVector:
    carry_counter = 0

    def __init__(self, entries):
        self.entries = entries

    def __xor__(self, other):
        print("WARNIGN")
        assert len(self.entries) == len(other.entries)
        result = []
        for i in range(len(self.entries)):
            result.append((self.entries[i] + other.entries[i]) % 2)
        return FormalVector(result)

    # We must always enforce the result of an addition is a bitvector
    def __add__(self, other):
        assert len(self.entries) == len(other.entries)
        result = []
        for i in range(len(self.entries)):
            result.append(self.entries[i] + other.entries[i])
        return FormalVector(result)

    def __rshift__(self, other):
        result = []
        for i in range(len(self.entries) - other):
            result.append(copy.copy(self.entries[i + other]))
        for i in range(other):
            result.append(LinearCombination([]))
        return FormalVector(result)

    def __lshift__(self, other):
        result = []
        for i in range(other):
            result.append(LinearCombination([]))
        for i in range(len(self.entries) - other):
            result.append(copy.copy(self.entries[i]))
        return FormalVector(result)

    def __mul__(self, other):
        result = make_const_vector(len(self.entries), 0)
        shift = 0
        while other != 0:
            if other & 1 == 1:
                result += (self << shift)
            shift += 1
            other >>= 1
        return result

    def rot_left(self, amount):
        # print((self << amount) + (self >> (len(self.entries) - amount)))
        # print()
        return (self << amount) + (self >> (len(self.entries) - amount))

    def with_mod_two_vars(self):
        carrys = make_var_vector(len(self.entries), f"c_{FormalVector.carry_counter}_", 'i')
        FormalVector.carry_counter += 1
        result = []
        for i in range(len(self.entries)):
            result.append(self.entries[i] + carrys.entries[i].scale(2))
        return FormalVector(result)

    def with_integer_carry_vars(self):
        carrys = make_var_vector(len(self.entries), f"c_{FormalVector.carry_counter}_", 'i')
        FormalVector.carry_counter += 1
        result = []
        for i in range(len(self.entries)):
            if i > 0:
                result.append(self.entries[i] + carrys.entries[i].scale(2) + carrys.entries[i - 1].scale(-1))
            else:
                result.append(self.entries[i] + carrys.entries[i].scale(2))
        return FormalVector(result)

    def get_bitvector_constraints(self):
        constraint_strings = []
        for m in self.entries:
            if len(str(m)) > 0:
                constraint_strings.append(str(m) + " >= 0")
                constraint_strings.append(str(m) + " <= 1")
        return constraint_strings

    def get_zero_and_bitvector_constraints(self, start, end):
        constraint_strings = []
        for i in range(0, len(self.entries)):
            m = self.entries[i]
            if len(str(m)) > 0:
                if start <= i < end:
                    constraint_strings.append(str(m) + " = 0")
                else:
                    constraint_strings.append(str(m) + " >= 0")
                    constraint_strings.append(str(m) + " <= 1")
        return constraint_strings

    def __str__(self):
        s = ""
        for comb in self.entries:
            s += str(comb) + "\n"
        return s


def make_var_vector(length, name, var_type):
    entries = []
    for i in range(length):
        v = Variable(name + str(i), var_type)
        entries.append(LinearCombination([Monomial(v, 1)]))
    return FormalVector(entries)


def make_const_vector(length, bits):
    while bits < 0:
        print("This shit is running")
        bits += (1 << length)
    entries = []
    for i in range(length):
        if bits & 1 == 1:
            entries.append(LinearCombination([Monomial(Variable.const, 1)]))
        else:
            entries.append(LinearCombination([]))
        bits >>= 1
    return FormalVector(entries)


def print_test_program():
    l = make_var_vector(64, 's', 'b')
    m = make_const_vector(64, 18446744071950099646)
    i1 = (l + m).with_integer_carry_vars()
    constraints = []
    constraints += i1.get_bitvector_constraints()
    out = (i1.rot_left(17) + l).with_integer_carry_vars()
    constraints += out.get_zero_and_bitvector_constraints(0, 54)

    print("Subject To")
    index = 0
    for c in constraints:
        print(f"c{index}: {c}")
        index += 1
    print(f"c{index}: ONE = 1")

    print("Bounds")
    for v, vtype in Variable.created_variables.items():
        if vtype == 'b':
            print(f"0 <= {v} <= 1")
        else:
            print(f"{v} free")

    print("Generals")
    for v, vtype in Variable.created_variables.items():
        print(v)

    print("End")


# (l, m) or (seedLo seedHi)
def xoroshiro(seed_tuple):
    l = seed_tuple[0]
    m = seed_tuple[1]
    s1 = m ^ l
    s0 = (l.rot_left(49) ^ s1) ^ (s1 << 21)
    s1 = s1.rot_left(28)
    return s0, s1


def write_program(constraints, path):
    f = open(path, "w")
    f.write("Subject To\n")
    index = 0
    for c in constraints:
        f.write(f"c{index}: {c}\n")
        index += 1
    f.write(f"c{index}: ONE = 1\n")

    f.write("Bounds\n")
    for v, vtype in Variable.created_variables.items():
        if vtype == 'b':
            f.write(f"0 <= {v} <= 1\n")
        else:
            f.write(f"{v} free\n")

    f.write("Generals\n")
    for v, vtype in Variable.created_variables.items():
        f.write(v + "\n")

    f.write("End\n")
    f.close()


def print_skull_program(skulls = 2):
    l = make_var_vector(64, 's', 'b')
    m = make_const_vector(64, 18446744071950099646)
    constraints = []
    skull_seed = xoroshiro(xoroshiro((l, m)))

    for i in range(skulls):
        print(i)
        l = skull_seed[0].with_mod_two_vars()
        m = skull_seed[1].with_mod_two_vars()
        constraints += l.get_bitvector_constraints()
        constraints += m.get_bitvector_constraints()
        i1 = (l + m).with_integer_carry_vars()
        constraints += i1.get_bitvector_constraints()
        out = (i1.rot_left(17) + l).with_integer_carry_vars()
        constraints += out.get_zero_and_bitvector_constraints(64-6, 64)
        skull_seed = xoroshiro(xoroshiro(xoroshiro((l, m))))

    write_program(constraints, f"{skulls}_wither_skulls_one_in_64.lp")


# public static long mixStafford13(long seed) {
#       seed = (seed ^ seed >>> 30) * -4658895280553007687L;
#       seed = (seed ^ seed >>> 27) * -7723592293110705685L;
#       return seed;
#   }
#Output not enforced as bit vector
def mix_stafford(constraints, seed):
    i1 = (seed ^ (seed >> 30)).with_mod_two_vars()
    constraints += i1.get_bitvector_constraints()
    i2 = (i1 * ((1 << 64)-4658895280553007687)).with_integer_carry_vars()
    constraints += i2.get_bitvector_constraints()
    i3 = (i2 ^ (i2 >> 30)).with_mod_two_vars()
    constraints += i3.get_bitvector_constraints()
    i4 = (i3 * ((1 << 64)-7723592293110705685)).with_integer_carry_vars()
    constraints += i4.get_bitvector_constraints()
    return i4 ^ (i4 >> 31)


def create_repeating_decorator_program():
    seed = make_var_vector(64, 's', 'b')
    offset = make_const_vector(64, (1 << 64)-7046029254386353131)
    i1 = (seed + offset).with_integer_carry_vars()
    constraints = i1.get_bitvector_constraints()

    l = mix_stafford(constraints, seed)
    m = mix_stafford(constraints, i1)

    constraints += l.get_bitvector_constraints()
    constraints += m.get_bitvector_constraints()
    i1 = (l + m).with_integer_carry_vars()
    constraints += i1.get_bitvector_constraints()
    out = (i1.rot_left(17) + l).with_integer_carry_vars()
    constraints += out.get_zero_and_bitvector_constraints(1, 60)

    x = xoroshiro((l, m))
    l = x[0].with_mod_two_vars()
    m = x[1].with_mod_two_vars()

    print(len(constraints))
    constraints += l.get_bitvector_constraints()
    constraints += m.get_bitvector_constraints()
    i1 = (l + m).with_integer_carry_vars()
    constraints += i1.get_bitvector_constraints()
    out = (i1.rot_left(17) + l).with_integer_carry_vars()
    constraints += out.get_zero_and_bitvector_constraints(1, 60)

    write_program(constraints, "decorator_double_repeat.lp")


def solve(path):
    m = gp.read(path)
    m.Params.PoolSolutions = 1
    m.Params.PoolSearchMode = 1
    m.Params.MIPFocus = 1

    m.optimize()
    seeds = []
    checkRes = []
    for solNum in range(m.getAttr("SolCount")):
        seed = [None] * 64
        m.Params.SolutionNumber = solNum
        for var in m.getVars():
            if var.VarName[0] == 's':
                seed[int(var.VarName[1:])] = str(int(var.Xn))
        seeds.append("".join(seed[::-1]))
    print(seeds)
    print([int(s, 2) for s in seeds])
    print("Solutions valid?: ", False not in checkRes)


if __name__ == '__main__':
    print_skull_program(20)
    path = ""
    solve(path)
    # create_repeating_decorator_program()
    # solve(r"C:\Users\Matthew\Documents\PythonWorkspace\xoroshiroilp\decorator_double_repeat.lp")
    # print_skull_program(12)
    #
    # low_seed = int("1000101010010000000100000011100001100000110000101110101011100000", 2)
    # intermediate = (low_seed + 18446744071950099646) % (2 ** 64)
    # intermediate = (intermediate >> (64 - 17)) | ((intermediate << 17) % (2 ** 64))
    # print(bin((intermediate + low_seed) % (1 << 64)))
