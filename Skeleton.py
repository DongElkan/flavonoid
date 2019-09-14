"""
Basic skeleton object
"""


class _Skeleton():
    """
    Initialize the skeleton project
    """
    __slots__ = ("ringA",
                 "ringB",
                 "ringC",
                 "name",
                 "siderings",
                 "corders",
                 "sideatoms",
                 "derivatives",
                 "label",
                 "connects")

    def __init__(self):
        """
        set defaulted None to each attribute defined
        in slots
        """
        for key in self.__slots__:
            setattr(self, key, None)

    def __setattr__(self, key, value):
        """
        set values to keys
        """
        setattr(self, key, value)
