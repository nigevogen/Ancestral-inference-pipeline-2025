def expand_path(path):
    """Converts a collapsed path list [AAA, AAT, TAT] to 
    an expanded list [(AAA, AAT), (AAT, TAT)] of
    codon tuples representing one-hit changes.
    
    Parameters
    ----------
    path : list of str
        Collapsed representation of the codon transition path, (AAA, AAT, TAT).
        The transition path is represented by start and end codons and the
        list of one-hit codons needed to reach the end codon.
    
    Returns
    -------
    list of tuple
        The expanded path is a list of tuples representing edges in the codon
        transition path with a one-hit change.

    """
    return list(zip(path[:-1], path[1:]))

def collapse_path(path):
    """Converts an extended path list [(AAA, AAT), (AAT, TAT)] to a
    collapsed path list [AAA, AAT, TAT].
    
    Parameters
    ----------
    path : list of tuple
        Extended representation of the codon transition path, 
        [(AAA, AAT), (AAT, TAT)].
    
    Returns
    -------
    list of str
        The collapsed path is a list of str, each str representing
        one-hit change nodes in the codon transition path.

    """
    return [p[0] for p in path[:-1]] + list(path[-1])
