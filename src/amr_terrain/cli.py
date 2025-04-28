from .backendInterface import amrBackend


def main():
    from sys import argv

    amrRef = amrBackend(argv[1])
    amrRef.createDomain()
    amrRef.createAMRFiles()
