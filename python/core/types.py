import libjad_DelphesAnalysis as j

def StringVector(iterable):
    i=j.jad_StringVector()
    i.extend(list(iterable))
    return i

def DoubleVector(iterable):
    i=j.jad_DoubleVector()
    i.extend(list(iterable))
    return i

def IntVector(iterable):
    i=j.jad_IntVector()
    i.extend(list(iterable))
    return i

def FileObjectVector(iterable):
    i=j.jad_FileObject()
    i.extend(list(iterable))
    return i

def FileMapVector(iterable):
    i=j.jad_FileMap()
    i.extend(list(iterable))
    return i
