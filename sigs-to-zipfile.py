#! /usr/bin/env python
# modified from: https://github.com/dib-lab/sourmash/pull/1349#issuecomment-787482809
import sys
import zipfile
import sourmash
import argparse

def main():
    p = argparse.ArgumentParser()
    p.add_argument('zipfile')
    p.add_argument('signatures', nargs='*')
    p.add_argument('--sig-pathlist')
    args = p.parse_args()

    zf = zipfile.ZipFile(args.zipfile, 'w')

    siglist = [x.rstrip() for x in open(args.sig_pathlist)]
    all_sigs = siglist + args.signatures

    n = 0
    for filename in all_sigs:
        print(f"reading signatures from '{filename}'")
        for sig in sourmash.load_file_as_signatures(filename):
            # I think the purpose of this is that zip needs something unique for each file.
            # Can't (generally) trust filenames, so use md5sum instead?
            md5 = 'signatures/' + sig.md5sum() + '.sig'
            sigstr = sourmash.save_signatures([sig])
            zf.writestr(md5, sigstr)
            n += 1

    print(f"wrote {n} signatures to '{args.zipfile}'")

    return 0


if __name__ == '__main__':
    sys.exit(main())
