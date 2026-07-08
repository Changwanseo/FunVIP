#!/usr/bin/env python
"""Generate docs/parameters.md from the live argparse parser in command.py.

Run from the repo root:  python docs/gen_param_reference.py > docs/parameters.md
This keeps the parameter reference in sync with funvip/src/command.py.
"""
import sys

sys.argv = ["FunVIP"]  # so the parser's parse_args() call sees no arguments

from funvip.src.command import CommandParser


def _type_name(action):
    if action.const is True or action.nargs == 0:
        return "flag"
    if action.type is not None:
        return getattr(action.type, "__name__", str(action.type))
    if action.nargs in ("*", "+"):
        return "list"
    return "str"


def main():
    cp = CommandParser()
    cp.get_args()  # builds cp.parser (harmlessly parses no arguments)
    parser = cp.parser

    out = [
        "# FunVIP command-line parameters",
        "",
        "_Auto-generated from `funvip/src/command.py` by `docs/gen_param_reference.py`;",
        "do not edit by hand. Regenerate after changing the CLI._",
    ]
    for group in parser._action_groups:
        actions = [a for a in group._group_actions if a.option_strings]
        if not actions:
            continue
        title = (group.title or "options").strip()
        out += ["", f"## {title}", "", "| Option | Type | Description |", "|---|---|---|"]
        for a in actions:
            opts = ", ".join(f"`{o}`" for o in a.option_strings)
            help_ = (a.help or "").replace("|", r"\|").replace("\n", " ")
            out.append(f"| {opts} | {_type_name(a)} | {help_} |")
    print("\n".join(out))


if __name__ == "__main__":
    main()
