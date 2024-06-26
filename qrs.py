"""
Program that finds the smalles solution to a cover problem with a specific even-odd rule
for combining different sets.
"""

from typing import Iterator, List, Set


def choose_n_filtered_recurse(
    problem: List[List[bool]],
    column: int,
    n: int,
    index: int,
    selection: List[int],
    result: List[List[int]],
):
    """Recursive implementation for choose_n_filtered"""
    if n == 0:
        result.append(selection.copy())
        return

    if index == len(problem):
        return

    choose_n_filtered_recurse(problem, column, n, index + 1, selection.copy(), result)

    if problem[index][column]:
        selection.append(index)
        choose_n_filtered_recurse(
            problem, column, n - 1, index + 1, selection.copy(), result
        )


def choose_n_filtered(
    problem: List[List[bool]], column: int, n: int
) -> List[List[int]]:
    """
    Returns all ways of choosing N from a column of the problem, filtered by colum value

    Arguments
    ---------
    problem :
        The problem instance
    column :
        Which column to select from
    n :
        The number of items to select
    """
    result: List[List[int]] = []
    choose_n_filtered_recurse(problem, column, n, 0, [], result)
    return result


def find_mult_covers(problem: List[List[bool]], columns: List[int]) -> List[Set[int]]:
    """Find all row multiplications that covers a set of columns

    Arguments
    ---------
    problem :
        The problem instance
    columns :
        Which columns must be covered
    """

    selection: List[Set[int]] = []

    for column_idx, column in enumerate(columns):
        num_true = 0
        for row_idx in range(len(problem)):
            if problem[row_idx][column]:
                num_true += 1

        new_selection: List[Set[int]] = []

        # TODO: Store failures to fail sooner
        # TODO: Check if a set is already computed
        # TODO: The result should really be Set[Set[int]]

        for num in range(1, num_true + 1, 2):
            rows = choose_n_filtered(problem, column, num)

            if column_idx == 0:
                for row in rows:
                    new_selection.append(set(row))
                continue

            for sel in selection:
                for row in rows:

                    comb = sel.union(row)

                    count = 0
                    for r in comb:
                        if problem[r][columns[column_idx]]:
                            count += 1

                    if count % 2 == 0:
                        continue

                    if sel.issuperset(row):
                        new_selection.append(sel)
                        continue

                    are_compatible = True
                    for col in columns[: column_idx + 1]:
                        count = 0
                        for r in comb:
                            if problem[r][col]:
                                count += 1

                        if count % 2 == 0:
                            are_compatible = False
                            break

                    if are_compatible:
                        new_selection.append(comb)

        if len(new_selection) == 0:
            return []

        selection = new_selection

    return selection


def choose_n_recurse(
    options: List[int],
    n: int,
    index: int,
    selection: List[int],
    result: List[List[int]],
):
    """Recursive implementation of choose_n"""
    if n == 0:
        result.append(selection.copy())
        return

    if len(options) == index + n:
        for opt in options[index:]:
            selection.append(opt)

        result.append(selection.copy())
        return

    choose_n_recurse(options, n, index + 1, selection.copy(), result)

    selection.append(options[index])
    choose_n_recurse(options, n - 1, index + 1, selection.copy(), result)


def choose_n(options: List[int], n: int) -> List[List[int]]:
    """
    Returns all ways of choosing N from a collection of integers

    Arguments
    ---------
    options :
        The integers to choose from
    n :
        The number of items to select

    """
    result: List[List[int]] = []
    choose_n_recurse(options, n, 0, [], result)
    return result


def all_splits_recursive(
    count: int, splits: int, upper: int, current: List[int], result: List[List[int]]
):
    """Recursive implementation of all_splits"""
    if splits == 0:
        current.append(count)
        result.append(current.copy())
        return

    lower = count // (splits + 1)
    if count % (splits + 1) == 0:
        lower -= 1

    for i in range(min(upper, count - splits), lower, -1):
        new = current.copy()
        new.append(i)
        all_splits_recursive(count - i, splits - 1, i, new, result)


def all_splits(count: int, splits: int) -> List[List[int]]:
    """
    Returns all unique ways of splitting the number line up to a value into parts

    Arguments
    ---------
    count :
        The length of the number line
    splits :
        The number of splits to do (one less than the number of partitions)
    """
    result: List[List[int]] = []
    all_splits_recursive(count, splits, count, [], result)
    return result


def partitions_recursive(
    remaining: Set[int], splits: List[int], split_idx: int, current: List[List[int]]
) -> Iterator[List[List[int]]]:
    """"""
    if split_idx == len(splits):
        yield current
        return

    for choice in choose_n(list(remaining), splits[split_idx]):
        new_current = current.copy()
        new_current.append(choice)

        yield from partitions_recursive(
            remaining.difference(choice), splits, split_idx + 1, new_current
        )


def all_partitions(columns: int, partitions: int) -> Iterator[List[List[int]]]:
    """
    Returns all ways to partition the number of the number line into partitions

    Arguments
    ---------
    columns :
        The number of columns in the problem, i.e. the length of the number line
    partitions :
        The number of partitions to divide them into, i.e. number of splits plus one
    """
    for split in all_splits(columns, partitions - 1):
        yield from partitions_recursive(set(range(columns)), split, 0, [])


def find_smallest_covers(
    problem: List[List[bool]], find_all_solutions: bool = False
) -> List[List[List[Set[int]]]]:
    """Finds the smallest cover of the rows including multiplications.

    Arguments
    ---------
    problem :
        The problem instance
    find_all_solutions :
        Whether to print the first solution or all solutions with the same fewest number
        of multiplications
    """
    columns = len(problem[0])

    for i in range(columns):
        all_solutions: List[List[List[Set[int]]]] = []
        for partition in all_partitions(columns, i + 1):

            selections: List[List[Set[int]]] = []
            for cols in partition:
                cover = find_mult_covers(problem, cols)
                if len(cover) == 0:
                    break

                selections.append(cover)

            if len(selections) == len(partition):
                if find_all_solutions:
                    all_solutions.append(selections)
                else:
                    return [selections]

        if len(all_solutions) > 0:
            return all_solutions

    return []


problem = [
    [False, True, True, False, True],
    [False, False, True, True, False],
    [False, True, False, False, False],
    [True, True, False, False, False],
]
# problem = [
#     [True, False, True, True],
#     [False, False, True, False],
#     [False, True, True, False],
#     [False, False, True, True],
# ]

# problem = [
#     [False, True, True, False],
#     [True, True, False, True],
# ]

# problem = [
#     [False, True, True, False, True, False, True],
#     [False, True, False, True, True, True, False],
#     [True, True, False, True, False, False, True],
# ]
#
# problem = [
#     [False, True, True, True, False, False, True],
#     [True, False, True, True, False, True, False],
#     [False, True, True, False, True, True, False],
# ]

solution = find_smallest_covers(problem, find_all_solutions=False)
print(solution)
