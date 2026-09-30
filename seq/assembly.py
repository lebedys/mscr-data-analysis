"""The a-b-c-a-b-c assembly: fitting called parts onto the slot scheme."""

# the cassette repeats, so the slot after c opens the next one
NEXT_SLOT = {"a": "b", "b": "c", "c": "a"}

RESIDUAL = "residual"


def decompose(parts, slots):
    """First a-b-c-a-b-c assignment as (slot, option) pairs; None if the list does not fit."""

    def consume(index, pattern, assigned):
        if not pattern:
            return assigned if index == len(parts) else None
        for option in slots[pattern[0]]:
            end = index + len(option)
            if parts[index:end] == option:
                found = consume(end, pattern[1:], assigned + [(pattern[0], tuple(option))])
                if found:
                    return found
        return None

    return consume(0, "abcabc", [])


def assign_slots(parts, slots):
    """One slot per part occurrence, RESIDUAL where no legal chain covers it.

    Accepts any number of cassettes and unmatched flanking parts; maximises parts covered.
    """

    def covered(run):
        return sum(end - start for _, start, end in run)

    best = []
    for first_index in range(len(parts)):
        for first_slot in NEXT_SLOT:
            stack = [(first_index, first_slot, [])]
            while stack:
                index, slot, assigned = stack.pop()
                if covered(assigned) > covered(best):
                    best = assigned
                for option in slots[slot]:
                    end = index + len(option)
                    if parts[index:end] == option:
                        stack.append((end, NEXT_SLOT[slot], assigned + [(slot, index, end)]))

    placed = {index: slot for slot, start, end in best for index in range(start, end)}
    return [placed.get(index, RESIDUAL) for index in range(len(parts))]
