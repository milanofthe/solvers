"""Check every DOI in the library against Crossref.

A reference database is only worth the references. A DOI that resolves to a
different paper than the one printed beside it is the kind of error nobody
notices, because it looks right and nobody clicks it, so it is checked here
instead: the DOI is resolved, and the title Crossref returns is compared with
the title in the file.

The comparison is loose on purpose. Publishers hold titles with different
capitalisation, different dashes and the occasional subtitle, so what is
compared is the sequence of alphanumeric words. What matters is catching a DOI
that points at another paper, not one that points at the same paper spelled
differently.

Books and reports are often not registered with Crossref at all, and a DOI that
cannot be resolved is reported separately from one that resolves to the wrong
thing.

Usage:
    python tools/check_references.py
"""

import io
import json
import os
import re
import sys
import time
import urllib.error
import urllib.request

HERE = os.path.dirname(os.path.abspath(__file__))
METHODS = os.path.join(HERE, "..", "methods")
AGENT = "solvers-reference-check (https://github.com/milanofthe/solvers)"
WORD = re.compile(r"[a-z0-9]+")


def words(title):
    return WORD.findall(title.lower())


def overlap(left, right):
    """How much of the shorter title the longer one contains."""
    a, b = set(words(left)), set(words(right))
    if not a or not b:
        return 0.0
    return len(a & b) / min(len(a), len(b))


def crossref(doi):
    request = urllib.request.Request(
        f"https://api.crossref.org/works/{urllib.parse.quote(doi)}",
        headers={"User-Agent": AGENT},
    )
    with urllib.request.urlopen(request, timeout=30) as response:
        payload = json.load(response)
    titles = payload["message"].get("title") or []
    return titles[0] if titles else ""


def references():
    """Every distinct DOI in the library, with the title it is filed under and
    the methods that cite it."""
    found = {}
    for folder in sorted(os.listdir(METHODS)):
        directory = os.path.join(METHODS, folder)
        if not os.path.isdir(directory):
            continue
        for name in sorted(os.listdir(directory)):
            entry = json.load(io.open(os.path.join(directory, name), encoding="utf-8"))
            for reference in entry.get("references", []):
                doi = reference.get("doi")
                if not doi:
                    continue
                record = found.setdefault(doi, {"title": reference.get("title", ""), "ids": []})
                record["ids"].append(entry["id"])
    return found


def main():
    found = references()
    print(f"{len(found)} distinct DOIs")
    unresolved, mismatched = [], []
    for index, (doi, record) in enumerate(sorted(found.items()), start=1):
        try:
            registered = crossref(doi)
        except urllib.error.HTTPError as error:
            unresolved.append((doi, record, f"HTTP {error.code}"))
            continue
        except Exception as error:  # network, timeout, malformed payload
            unresolved.append((doi, record, type(error).__name__))
            continue
        score = overlap(record["title"], registered)
        if score < 0.6:
            mismatched.append((doi, record, registered, score))
        print(f"  [{index}/{len(found)}] {doi} {score:.2f}")
        time.sleep(0.2)

    if mismatched:
        print("\nDOIs that resolve to something else:")
        for doi, record, registered, score in mismatched:
            print(f"  {doi}  ({score:.2f} of the words match)")
            print(f"    file:     {record['title']}")
            print(f"    crossref: {registered}")
            print(f"    cited by: {', '.join(record['ids'][:6])}")
    if unresolved:
        print("\nDOIs Crossref does not hold (books and reports usually):")
        for doi, record, reason in unresolved:
            print(f"  {doi}  {reason}  {record['title'][:70]}")
    if not mismatched and not unresolved:
        print("\nevery DOI resolves to the paper it is filed under")
    return 1 if mismatched else 0


if __name__ == "__main__":
    import urllib.parse

    sys.exit(main())
