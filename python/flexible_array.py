from collections import OrderedDict

class FlexibleArray:
    def __init__(self):
        self.data = OrderedDict()
        self.next_key = 0  # Tracks next available key

    def append(self, obj):
        """Append an object to the next available integer key."""
        self.data[self.next_key] = obj
        self.next_key += 1

    def drop_first(self):
        """Remove 'count' oldest elements (prefix)."""
        #print(f"drop first")
        #print(self.data.keys())
        self.data.pop(next(iter(self.data)))  # Remove first item
        #print(self.data.keys())
        #print("----")

    def drop_last(self):
        """Remove 'count' newest elements (suffix)."""
        #print(f"drop last")
        #print(self.data.keys())
        #self.next_key = self.next_key - count
        #for _ in range(min(count, len(self.data))):
        self.data.popitem()  # Remove last item
        self.next_key = self.next_key - 1
        #print(self.data.keys())
        #print("----")

    def first(self):
        """Access the first element without knowing its key."""
        if self.data:
            return next(iter(self.data.values()))  # First element
        return None

    def last(self):
        """Access the last element without knowing its key."""
        if self.data:
            return next(reversed(self.data.values()))  # Last element
        return None

    def size(self):
        """Return the number of elements currently in the flexible array."""
        return len(self.data)

    def first_index(self):
        """Return the index of the first element, or None if empty."""
        if self.data:
            return next(iter(self.data.keys()))
        return None

    def last_index(self):
        """Return the index of the last element, or None if empty."""
        if self.data:
            return next(reversed(self.data.keys()))
        return None

    def __repr__(self):
        return repr(self.data)

# Example Usage:
#fa = FlexibleArray()
#fa.append("A")
#fa.append("B")
#fa.append("C")
#
#print(fa.first())        # Output: "A"
#print(fa.last())         # Output: "C"
#print(fa.size())         # Output: 3
#print(fa.first_index())  # Output: 0
#print(fa.last_index())   # Output: 2
#
#fa.drop_first(1)  # Remove "A"
#print(fa.first_index())  # Output: 1
#
#fa.drop_last(1)  # Remove "C"
#print(fa.last_index())   # Output: 1
