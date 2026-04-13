import os
import re

from setuptools import setup

LINK_PATTERN = r"(\[.*?\]\()([^http][^)]*)\)"
IMAGE_SRC_PATTERN = r'(<img\s+[^>]*src=")([^"]*)(")'
BASE_URL = "https://github.com/bioinform/somaticseq"


def _read_version() -> str:
    # Read __version__ from the _version.py file.
    version_file = os.path.join("somaticseq", "_version.py")
    namespace: dict[str, str] = {}
    with open(version_file, encoding="utf-8") as f:
        exec(f.read(), namespace)
    return namespace["__version__"]


def modify_markdown_for_3rd_party(base_markdown: str, version: str) -> str:
    def _replace_link(match: re.Match) -> str:
        # Replace relative links in .md from [text](RELATIVE/LINK) into
        # [text]({BASE_URL}/blob/TAG/RELATIVE/LINK)
        text = match.group(1)
        url = match.group(2)
        return f"{text}{BASE_URL}/blob/v{version}/{url})"

    def _replace_src(match: re.Match) -> str:
        # Replace relative image links.
        prefix = match.group(1)
        url = match.group(2)
        suffix = match.group(3)
        return f"{prefix}{BASE_URL}/raw/v{version}/{url}{suffix}"

    with_abs_url = re.sub(LINK_PATTERN, _replace_link, base_markdown)
    with_abs_img_src = re.sub(IMAGE_SRC_PATTERN, _replace_src, with_abs_url)
    return with_abs_img_src


with open("README.md", encoding="utf-8") as f:
    version = _read_version()
    description_for_3rd_party = modify_markdown_for_3rd_party(f.read(), version)


setup(
    long_description=description_for_3rd_party,
    long_description_content_type="text/markdown",
)
