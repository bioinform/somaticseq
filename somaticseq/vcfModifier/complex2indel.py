#!/usr/bin/env python3

def translate(refbase, altbase):
    
    offset  = 0
    
    if len(refbase) == len(altbase):
        return False
    
    elif len(refbase) == 1 or len(altbase) == 1:
        return ( (refbase, altbase), offset )
        
    else:
        for base_i, base_j in zip(refbase[::-1], altbase[::-1]):
            if base_i == base_j and (len(refbase) >= 2 and len(altbase) >= 2):
                refbase = refbase[:-1]
                altbase = altbase[:-1]
            else:
                break

        for base_i, base_j in zip(refbase, altbase):
            if base_i == base_j and (len(refbase) >= 2 and len(altbase) >= 2):
                refbase = refbase[1:]
                altbase = altbase[1:]
                offset += 1
            else:
                break

        return ( (refbase, altbase), offset )


if __name__ == '__main__':
    translate(refbase, altbase)
