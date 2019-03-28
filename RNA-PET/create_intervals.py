

def create_interval(pos, pos_type, dist):
    """ Creates a zero-based interval around the provided position (must also be 
        zero-based) of size dist on either side. """

    if pos_type == "start":
        start = pos - dist
        end = pos + dist + 1
    elif pos_type == "end":
        start = pos - dist - 1
        end = pos + dist

    return start, end

if __name__ == '__main__':
    main()
