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
    p.add_argument('--compression', type=int, default=9)
    p.add_argument('--ksize', type=int) #nargs="*", default=[21,31,51]) #allow multiple??
    args = p.parse_args()

    zf = zipfile.ZipFile(args.zipfile, 'w')

    siglist = [x.rstrip() for x in open(args.sig_pathlist)]
    all_sigs = siglist + args.signatures

    n = 0
    all_md5=set()
    for i, filename in enumerate(all_sigs):
        if n % 10000 == 0:
            print(f"... processing {n}th signature; currently reading signatures from '{filename}'")

        for sig in sourmash.load_file_as_signatures(filename, ksize=args.ksize):
            # zip needs a unique name for each signature. Use sig md5sum.
            md5= sig.md5sum()
            # if this is a duplicate md5sum, add _{number} to make it unique.
            if md5 in all_md5:
                sys.stderr.write(f"{str(sig)} has an md5sum identical to one already in the zipfile ({md5})")
                i=0
                full_md5 = f"{md5}_{i}"
                while full_md5 in all_md5:
                    i+= 1
                    full_md5 = f"{md5}_{i}"
                md5=full_md5
                sys.stderr.write(f"...adding unique md5 {md5} instead")

            all_md5.add(md5)
            md5_name = 'signatures/' + md5 + '.sig'
            sigstr = sourmash.save_signatures([sig], compression=args.compression)
            zf.writestr(md5_name, sigstr)
            n += 1

    print(f"wrote {n} signatures to '{args.zipfile}'")

    return 0


if __name__ == '__main__':
    sys.exit(main())
