## CADET Documentation

Install build dependencies:

```
pip install -r requirements.txt
```

Build the documentation locally from the `doc` folder:

```
myst build
```

The output is in `_build/html/` and can be opened with any browser.

Any changes to the documentation will automatically be pushed to the github-pages repository (<https://github.com/cadet/cadet.github.io>) using GitHub Actions.
