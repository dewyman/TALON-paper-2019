
def create_end_piece(pos, strand, dist):
    """ Creates a zero based interval of length 'dist' that starts inside the
        transcript and ends with the transcript end."""

    if strand == "+":
        interval_start = pos - dist
        interval_end = pos
    elif strand == "-":
        interval_start = pos
        interval_end = pos + dist
    else:
        raise ValueError("Strand must be '+' or '-'.")

    return interval_start, interval_end

def create_interval(pos, pos_type, dist):
    """ Creates a zero-based interval around the provided position (must also be 
        zero-based) of size dist on either side. """

    if pos_type == "left":
        start = pos - dist
        end = pos + dist + 1
    elif pos_type == "right":
        start = pos - dist - 1
        end = pos + dist

    return start, end

if __name__ == '__main__':
    main()
