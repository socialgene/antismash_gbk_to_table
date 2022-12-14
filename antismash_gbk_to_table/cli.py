#!/usr/bin/env python
import argparse
from pathlib import Path
from antismash_gbk_to_table.funcs import parse_and_write

parser = argparse.ArgumentParser(description="Parse antismash GenBank files")
parser.add_argument(
    "-i",
    "--input",
    metavar="filepath",
    help="Path to GenBank file",
    required=True,
)
parser.add_argument(
    "-o",
    "--output",
    metavar="filepath",
    help="File path to write to",
    required=True,
)
parser.add_argument(
    "-a",
    "--append",
    metavar="bool",
    default=False,
    help="numoutfiles file",
    required=False,
    action=argparse.BooleanOptionalAction,
)

parser.add_argument(
    "--header",
    metavar="bool",
    default=False,
    help="Write column headers?",
    required=False,
    action=argparse.BooleanOptionalAction,
)


def main():
    args = parser.parse_args()
    if Path(args.output).exists() and not args.append:
        raise FileExistsError(Path(args.output))
    parse_and_write(
        gbk_path_list=args.input,
        outpath=args.output,
        append=args.append,
        header=args.header,
    )


if __name__ == "__main__":
    main()

if __name__ == "__main__":
    sys.exit(main())  # pragma: no cover
