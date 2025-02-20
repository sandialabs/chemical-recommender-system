# Contributing Guidelines
### Fork CRS

*   If you have not already done so, create a fork of CRS on GitHub under your username.
    *   Sign on (via web) to https://github.com/sandialabs/chemical-recommender-system
    *   Make sure you are signed in to github
    *   Click on the 'Fork' button near the top right of the page.
*   Clone your fork of CRS with
    *   `git clone git@github.com:<username>/chemical-recommender-system`
    *   Or: `git clone https://github.com/<username>/chemical-recommender-system`
*   Each time you clone your fork,
    *   `git remote add upstream git@github.com:sandialabs/chemical-recommender-system` to add the original CRS repository as the `upstream` remote.
    *   Or: `git remote add upstream https://github.com/sandialabs/chemical-recommender-system`

### Update the Main Development Branch

To keep your `master` branch up-to-date with `upstream`:

*   `git fetch --all`
*   `git checkout master`
*   `git merge upstream/master`
*   `git push origin master`

You want to do this before starting work on a new feature branch.

### Create a Feature Branch

Create a local branch off of `master` on which to make your changes:

*   `git checkout master`
*   `git checkout -b <branchName>`

`<branchName>` can be whatever you like, though we have some recommendations:
*   Make the branch name descriptive; that is, avoid `fixSomeStuff`, `performanceTweaks`, and generic names along those lines.
*   To indicate your branch is intended solely for your own use, preface the branch name with your username, as in `<username>/<restOfBranchName>`.

### Make Your Changes

Do whatever work is necessary to address the issue you're tackling,
breaking your work into logical, compilable commits.  Feel free to
commit small chunks early and often in your local repository and then
use `git rebase -i` to reorganize your commits before sharing.  Make
sure the commit messages you will be sharing reference the appropriate
GitHub issue numbers.

### Update Your Branch

While working on your feature in your local `<branchName>` branch,
other commits will likely make it into the real CRS `master`
branch.  There are a variety of ways to merge these changes into your
local feature branch.  One possibility is

*   `git checkout <branchName>`
*   `git fetch --all`
*   `git merge upstream/master`

though there are others that are equally valid.

### Create a Pull Request

When your changes are ready to be integrated into CRS' `master` branch:

*   Push your local feature branch up to your fork with `git push -u origin <branchName>`.

*   Navigate to your fork of CRS on GitHub and create a new pull request:

*   Be sure you choose:
    *   base fork:  `sandialabs/chemical-recommender-system`
    *   base:  `master`
    *   head fork:  `<username>/chemical-recommender-system`
    *   compare:  `<branchName>`

### Feedback

At this point you'll enter into a stage where you and various CRS
developers will iterate back and forth until your changes are in an
acceptable state and can be merged in.  If you need to make changes to
your pull request, make additional commits on your `<branchName>`
branch and push them up to your fork.  Make sure you don't delete your
remote feature branch or your fork of CRS before your pull request
has been merged.

### Acknowledgement
Based on the `CONTRIBUTING.md` document from Trilinos.