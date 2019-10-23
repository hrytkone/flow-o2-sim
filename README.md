# ToyFlow

Some points to keep in mind:
- pT-dependence function does not correspond the real pT-dependence, it just looks like the vn(pT)
- REMEMBER to give percentage of lost particle in non-uniform case in range 0-1, not in 0-100!

# Running in puck
Useful commands to know for running in puck:
- Initializing the root: source setup.sh
- Check ejobs.sh for correct settings like how many jobs will be launched.
- Remember to compile the code if changes have been made.
- Launching massive jobs: ./ejobs.sh \<arguments\>
- Checking user job queue: squeue -u \<username\>
- Cancelling user jobs: scancel -u \<username\>
