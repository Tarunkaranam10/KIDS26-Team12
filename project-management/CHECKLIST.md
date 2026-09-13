# Team Lead Checklist

Use this page to get the team moving. It is intentionally short: a three-day project needs enough structure to coordinate work, not a second project to maintain.

## Start Here

1. **Access your assigned team repository.** Organizers will provide the repository and add team members. Confirm that you can open it on GitHub and clone it to your computer. Follow [Git and GitHub basics](../docs/git-github-basics.md) for the full path.
2. **Complete the project profile.** Agree on the question, inputs, expected output, tools, and team roles before pursuing a large implementation.
3. **Make a small first change.** Create a branch, update this README or document a data source, commit the change, and open a pull request. Use the terminal steps in [Git and GitHub basics](../docs/git-github-basics.md) or [GitHub Desktop](https://desktop.github.com/) if you prefer a graphical interface.
4. **Ask for help early.** Record blockers in an issue, raise them at a check-in, or ask a mentor. See [troubleshooting](../docs/troubleshooting.md) for common recovery steps.

## The Team

- **Members and roles:** Record these in [team.md](team.md).
- **Ways of working:** Use the [team lead checklist](CHECKLIST.md) to agree on branches, reviews, communication, and check-ins.
- **Current plan:** Keep small first tasks and risks in [project-plan.md](project-plan.md).

## Project Structure

Use the folders that fit your project. You do not need to fill every folder. The following is just a suggestion, yours might look different.

```text
data/raw/           Original inputs; do not edit in place
data/processed/     Cleaned or transformed data
models/             Models, predictions, or model notes
project-management/ Team plan, roles, decisions, and check-ins
src/                Reusable code, organized by purpose
docs/               Optional learning and troubleshooting guides
assets/             Images or other supporting project assets
```

## Data, meta data and secrets
**Do not commit passwords, API keys, private information, or identifiable human or clinical data. Check the source and license before sharing external data or media.**

## Resources

- New to Git or GitHub or need to know how to work with git in a shared repo: read [Git and GitHub basics](../docs/git-github-basics.md).
- Using Copilot agents: read [AI assistance](../docs/ai-guidance.md).
- Stuck during setup: open [troubleshooting](../docs/troubleshooting.md).
- Collaborating on changes: see [Contributing to your team](#contributing-to-your-team) below.

## Contributing to your team
### A Simple Workflow

1. Pick a small task or write down a blocker.
2. Create a branch with a clear name, such as `add-project-profile` or `fix-data-path`.
3. Make one focused change and commit it with a short message.
4. Push the branch and open a pull request.
5. Ask another teammate to look at the change before merging.
6. Update the README or project notes when the change affects how someone uses the project.

The [Git and GitHub basics](../docs/git-github-basics.md) guide explains each step, including a GitHub Desktop workflow.

### A Pull Request Is Ready When

- The change has a clear purpose.
- A teammate can understand what changed.
- You have recorded how you checked it, or explained why checking was not possible.
- Relevant assumptions, data sources, and limitations are documented.
- The change does not include credentials or sensitive data.

Small, incomplete pull requests are welcome when they make the current state visible and clearly describe what remains.


## Reproducibility and Attribution

Make work easier to inspect and reuse by keeping inputs, decisions, methods, and limitations visible. Prefer small readable steps over unexplained one-off commands. Cite data, code, models, and external resources that your project depends on. These practices help the next person understand what happened.

This repository is a reusable template. See [LICENSE.md](../LICENSE.md) for the licensing terms and update the project profile and attribution when you create a team project.

## Code of Conduct

In this repository, we use St. Jude's Code of Conduct document, outlining our expectations for all participants. Ask for help early, give feedback about the work rather than the person, and make room for different levels of experience.

For details, please visit: https://issuu.com/sjcrh/docs/st._jude_code_of_conduct.


## Team Leads: Before the Event

- [ ] Complete the [project profile](../README.md#project-profile).
- [ ] Agree on one communication channel and a short check-in rhythm.
- [ ] Create three to six small first tasks in the project board or [project-plan.md](project-plan.md).
- [ ] Use the plan and the expected output to suggest practical roles in [team.md](team.md).


## During the Three Days

### Day 1: Align and start

- Confirm the question, problem, or opportunity.
- Confirm the inputs and expected output.
- Make sure everyone can clone the repository and make a small change.
- Agree on branch, commit, and review habits.

### Day 2: Build and learn

- Keep tasks small enough to finish or review in one sitting.
- Record decisions that change the approach in a decision log (create `decisions.md` if useful).
- Document data sources, assumptions, and unexpected limitations as they appear.
- Check in briefly and redistribute work when someone is blocked.

### Day 3: Explain and hand off

- Decide what the final demo and booth must show as a team.
- Make the main workflow understandable to someone who was not in the room.
- Capture what worked, what did not, and what should happen next.
- Run the available checks and record their results.

# Final Output and Handoff

Use this space for the material that helps someone understand the project after the event.

- **Final demo or report:** The precomputed [Shiny app](../app/app.R) and final figures/results will be linked here during Day 3; current handoff status is in [verification status](../docs/14_VERIFICATION_STATUS.md).
- **Main result:** Pending biological fitting. The completed pre-event result is an auditable cohort/feature/model plan with a 7,707-candidate metadata cohort, 384,640-probe HM450/EPIC-v1 bridge, and leakage-safe elastic-net implementation.
- **How to reproduce or run it:** [Technical README](../README.md), [data acquisition guide](../docs/03_DATA_ACQUISITION.md), and [three-day runbook](../docs/05_THREE_DAY_HACKATHON_RUNBOOK.md).
- **Data and source notes:** [Acquisition provenance](../docs/03_DATA_ACQUISITION.md), [literature evidence](../docs/01_LITERATURE_EVIDENCE.md), and [feature-bridge provenance](../docs/17_450K_EPIC_FEATURE_BRIDGE.md).
- **Known limitations:** No biological model has yet been trained; the full historical matrix path/checksum and R runtime remain to be verified; canonical LOH/TAI cannot be calculated from ordinary methylation total CN; PBTP use requires the internal readiness gate.
- **Next steps:** Verify the reported full data download, run R smoke tests, extract the frozen feature matrix, fit the real-data baseline/null, and freeze before locked CNS/PBTP evaluation.

Keep generated figures and reports clearly named. Do not commit sensitive data or files that cannot be redistributed.

## Communications

Keep communication easy to find and easy to use during the three-day event.

- **Primary channel:** [Team 12 Slack](https://stjudebiohackathon.slack.com/archives/C0BSC28M3U6)
- **Slack team channel:** [Team 12 Slack](https://stjudebiohackathon.slack.com/archives/C0BSC28M3U6)
- **Team leads:** [Evan Savage (@esavage111)](https://github.com/esavage111) and [Susanna Downing (@sdowning12)](https://github.com/sdowning12)
- **Mentor or support contact:** Evan Savage (@esavage111)
- **Check-in time:** 09:00 planning, 12:00 milestone/cut review, and 16:30 wrap-up during the event; adjust in the team channel.
- **Slack general channel:** [St. Jude Biohackathon general Slack](https://stjudebiohackathon.slack.com/archives/C04JD4M3TCM)

Use `project-management/check-in.md` for short updates when useful (create the file if needed). Do not store private contact details or sensitive project information in this public repository.

## Optional Templates

- [Team and roles](team.md)
- [Project plan](project-plan.md)

Use only the templates that help. The repository should make progress easier, not require paperwork for its own sake.
