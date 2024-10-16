## Advanced usage and tips using FunVIP

### Preparing Database
Preparing database is quite tricky when the reference does not provides supplementary data formatted in table

#### If table is possible from web viewer
Try to use "From web" function of Microsoft Excel
![image](https://github.com/user-attachments/assets/f715a86f-4195-4622-a768-ab0cc753d610)

#### If only pdf is available
Use ABBYY Finereader pdf if you have enough funding
![image](https://github.com/user-attachments/assets/5f504692-88f6-41ac-ac16-8372591c8307)

An alternative open source option is camelot. It is quite tricky to install, but better than writing all values one by one

[How to use camelot](https://camelot-py.readthedocs.io/en/master/user/cli.html)


### Use FunVIP iteratively
Utilize validation ability of FunVIP as much as possible.

Try to put all putative database and queries that you are going to use for the analysis at the initial run.

Run initial run with fast mode (because there will be a lot of sequneces).

Filter out invalid databases according to the result of the first run.

Do this iteratively until your dataset seems to be fine.

At the last run, run with the accurate mode.

### SVG files are good at searching and editing
The visualized tree are exported in svg format. 

They can be edited with Microsoft Powerpoint, Adobe Illustrator/Indesign, Affinity Designer/Publisher, and Inkscape.

Also, using internet browser such as chrome, you can search trees with Ctrl + F


### Change the criteria using --continue
FunVIP --continue option enables running FunVIP again with other options!

1. Copy previous FunVIP run for backup
2. From the previous command, add this command
```
--continue --outdir <location_of_the_previous_output> --runname <runname_of_previous_output> --step <step_you_want_to_restart>
```
3. Run it again with different criteria!


