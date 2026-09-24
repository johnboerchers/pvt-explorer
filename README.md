# p-v-T Surface Explorer

![p-v-T explorer showing a 1 MPa isobar crossing the vapor dome of water](Images/screenshot_p.png)

An interactive 3D picture of how a pure substance (like water) behaves as you change its **pressure (p)**, **specific volume (v)**, and **temperature (T)**. You can rotate the surface, slice it at constant p, v, or T, and see how the familiar p-v, T-v, and p-T diagrams from class are just slices and "shadows" of this one 3D shape.

**No coding is required.** You only need to install Python once and type one command to run the program. This guide walks you through every step.

*Instructors: see [For instructors](#for-instructors) at the end for teaching notes, accuracy, and limitations.*

---

## Quick start (if you already have Python)

Open a terminal in the folder containing `pvt_explorer.py` and run:

```
python -m pip install numpy scipy matplotlib CoolProp
python pvt_explorer.py
```

If that worked, skip to [Using the program](#using-the-program). If not, read on.

---

## Step 1: Install Python

We recommend **Anaconda**, a free bundle that includes Python and most of the scientific packages this program needs.

1. Go to **https://www.anaconda.com/download** and download the installer for your computer (Windows or Mac).
2. Run the installer and accept the default options.
   - **Windows:** when asked, installing "Just Me" is fine.
3. When it finishes, you're done. You don't need to open anything yet.

> **Already have Python from python.org?** That works too. On Windows, make sure you checked the box **"Add Python to PATH"** during installation. If you didn't, reinstall and check it.

---

## Step 2: Open a terminal

A terminal is a window where you type commands instead of clicking.

- **Windows (with Anaconda):** press the Start button, type **Anaconda Prompt**, and open it.
- **Windows (python.org install):** press Start, type **cmd**, and open **Command Prompt**.
- **Mac:** press `Cmd + Space`, type **Terminal**, and press Enter.

You'll see a line ending in `>` or `$` with a blinking cursor. That's where you type.

---

## Step 3: Install the packages

The program uses four add-on packages:

| Package | What it does |
|---|---|
| `numpy` | Math with arrays of numbers |
| `scipy` | Equation solving |
| `matplotlib` | Plotting and the interactive window |
| `CoolProp` | Real property data for water and other fluids (same data as your steam tables) |

In the terminal, type this and press Enter:

```
python -m pip install numpy scipy matplotlib CoolProp
```

- **Mac:** if you get `command not found: python`, use `python3` instead of `python` everywhere in this guide.
- **Windows:** if you get `'python' is not recognized`, try `py` instead of `python`.

You'll see a lot of text scroll by. It's finished when you get your blinking cursor back and see a line starting with `Successfully installed` (or `Requirement already satisfied`, which just means you already had it).

> **If CoolProp won't install:** don't worry. The program still runs and automatically switches to a "van der Waals" model fluid, which has the same shape (vapor dome, critical point, etc.) but uses dimensionless units instead of real water data. Anaconda users can also try:
> ```
> conda install -c conda-forge coolprop
> ```

---

## Step 4: Download and find the program

1. On the GitHub page for this project, click the green **Code** button, then **Download ZIP**. Unzip it (on Windows: right-click → *Extract All*). You don't need a GitHub account.
2. Put the unzipped folder somewhere easy to find. We suggest renaming it `thermo` and putting it on your Desktop.
3. In the terminal, move into that folder using the `cd` ("change directory") command:

   **Windows:**
   ```
   cd Desktop\thermo
   ```
   **Mac:**
   ```
   cd Desktop/thermo
   ```

   **Tip:** type `cd ` (with a space after it), then drag the `thermo` folder from your file browser into the terminal window. It will fill in the full path for you. Press Enter.

4. Check that you're in the right place by listing the files:
   - Windows: `dir`
   - Mac: `ls`

   You should see `pvt_explorer.py` in the list. If instead you see another folder (unzipping sometimes creates a folder inside a folder), `cd` into that one too and check again.

---

## Step 5: Run it

```
python pvt_explorer.py
```

After a few seconds (it prints `Building surface ...` while it calculates), a window opens. **Maximize the window** so everything has room.

To close the program, close the window. You can run it again any time by repeating Steps 2, 4 (the `cd` part), and 5.

---

## Using the program

The window has three parts: the 3D surface on the left, three 2D diagrams on the right, and controls along the bottom.

### The 3D surface (left)
- **Click and drag** to rotate it. Look at it from different angles; try looking straight down each axis.
- Colors show the phase:
  - **Blue:** compressed liquid
  - **Green:** liquid + vapor mixture (inside the vapor dome)
  - **Orange:** superheated vapor / gas
  - **Gray:** supercritical fluid
- The **black line** is the saturation curve (the vapor dome), and the **black dot** is the critical point.
- The p and v axes are logarithmic (each tick is 10 times the previous one). Otherwise the liquid region would be too thin to see.

### The 2D diagrams (right)
These are the p-v, T-v, and p-T diagrams from your textbook. The thin gray lines are reference curves (isotherms, isobars, isochores).

### The controls (bottom)
- **Hold constant:** choose which property to keep fixed.
  - **p (isobar):** constant pressure, like boiling water in a pot or a piston-cylinder
  - **v (isochore):** constant specific volume, like heating a sealed rigid tank
  - **T (isotherm):** constant temperature
- **Slider:** drag it to change the value being held constant. The red plane in the 3D plot shows where you're slicing, and the **red curve** is the process line.
- The **red curve appears in all three 2D diagrams.** The diagram with the **red border** is the one the slice actually lives in; the other two show how that same process looks when projected onto other planes.
- **Wall projections:** shows dashed "shadows" of the dome and red curve on the back walls of the 3D box. Those shadows *are* the 2D diagrams.
- **Show surface:** hides the surface so you can see the dome and red curve clearly.
- **Text at the bottom:** describes what is physically happening along the red curve, with saturation values you can check against your steam tables.

---

## Things to try

1. **Hold T constant at 200 °C.** Read p_sat, v_f and v_g from the text at the bottom, then look them up in the saturated water table in your textbook. Do they match?
2. **Hold p constant** and slide it upward. What happens to the flat boiling segment in the T-v diagram as you approach the critical point? What happens above it?
3. **Hold T constant.** Why does the entire horizontal part of the isotherm in the p-v diagram become a single point in the p-T diagram?
4. **Hold v constant (rigid tank).** Try a value smaller than the critical specific volume, then one larger. Where does the tank end up: full of liquid or full of vapor?
5. Turn off **Show surface** and turn on **Wall projections**, then rotate the 3D plot until you're looking straight at one wall. What do you see?

---

## Optional extras

Other fluids (any fluid CoolProp supports):
```
python pvt_explorer.py --fluid R134a
python pvt_explorer.py --fluid CarbonDioxide
python pvt_explorer.py --fluid Nitrogen
```

Van der Waals model fluid (dimensionless "reduced" units):
```
python pvt_explorer.py --eos vdw
```

Save pictures instead of opening a window (creates `pics_p.png`, `pics_v.png`, `pics_T.png`):
```
python pvt_explorer.py --save pics
```

---

## Troubleshooting

**`'python' is not recognized` (Windows) or `command not found: python` (Mac)**
Try `py` (Windows) or `python3` (Mac). With Anaconda, make sure you opened **Anaconda Prompt**, not the regular Command Prompt.

**`No module named 'numpy'` (or scipy, matplotlib)**
The packages were installed for a different copy of Python than the one running the program. Run the install command again using exactly the same word (`python`, `py`, or `python3`) you use to run the program:
```
python -m pip install numpy scipy matplotlib CoolProp
```

**`can't open file 'pvt_explorer.py': No such file or directory`**
The terminal isn't in the folder with the program. Redo the `cd` step in Step 4 and check with `dir` / `ls`.

**It says `CoolProp not found ... using a van der Waals fluid instead`**
The program is still working, just with the model fluid. Try installing CoolProp again (see Step 3) if you want real water data.

**The window opens but the plot won't rotate**
Check the small toolbar at the bottom or top of the window. If the magnifying glass (zoom) or cross-arrows (pan) button is highlighted, click it again to turn it off.

**Using Spyder and the plot appears as a still image in the console**
Go to *Tools → Preferences → IPython console → Graphics*, set *Backend* to **Automatic**, click OK, and restart Spyder (or restart the kernel). Then run the file again.

**Using Jupyter Notebook or Google Colab**
Notebooks don't show the interactive window. Run the program from a terminal as described above. In Colab, you can still use `--save` to make pictures.

**The slider feels a little slow when dragging pressure**
That's normal. At each position the program calculates hundreds of property values. Drag slowly, or click a spot on the slider bar to jump to it.

**Things overlap or look cramped**
Maximize the window. The layout is designed for a full-size screen.

**Still stuck?**
Copy the *entire* error message (the last few lines are the most important) and bring it to office hours or email it to your instructor.

---

## For instructors

This tool was built for a sophomore-level engineering thermodynamics course. It is meant to help students connect the 3D p-v-T surface to the 2D property diagrams and tables they use every day. It is free to use, modify, and share under the MIT License.

### Topics it supports
- The vapor dome: saturated liquid and saturated vapor lines, the critical point, and the compressed liquid, saturated mixture, and superheated vapor regions
- Why pressure and temperature are not independent inside the dome, and why the two-phase region collapses to a single curve (the vaporization curve) in the p-T diagram
- Reading saturation properties (p_sat, T_sat, v_f, v_g) and checking them against steam tables
- Common processes: constant-pressure heating (piston-cylinder, boiler), constant-volume heating (rigid tank), and isothermal compression
- The rigid-tank question of whether a heated mixture ends as liquid or vapor, depending on whether v is below or above the critical specific volume
- Behavior above the critical point, and the approach to ideal-gas behavior at large specific volume
- (With `--eos vdw`) how a cubic equation of state produces a vapor dome, and the Maxwell equal-area construction

### Ways to use it
- **Live demo in lecture.** Project it and drag the slider while students predict what the red curve will do in each 2D diagram before you move it. Maximizing the window and turning off *Show surface* makes the dome easier to see from the back of the room.
- **In-class activity or homework.** The "Things to try" questions above work as a short worksheet. Asking students to verify values against the saturated water tables in the course textbook builds trust in both the tool and the tables.
- **Slides and exams.** `python pvt_explorer.py --save name` writes a PNG for each slicing mode without opening a window. To change the slice values used, edit the `defaults` in the `CoolPropFluid` class near the top of the script.
- **Other fluids.** Refrigerant cycles can use `--fluid R134a`. Any fluid name supported by CoolProp works.

### Accuracy
- Real-fluid properties come from [CoolProp](http://www.coolprop.org), which implements reference-quality Helmholtz-energy equations of state. For water this is the IAPWS-95 formulation, the same basis as modern steam tables, so values agree with textbook tables to the precision the tables print. Example: at 1 MPa the program gives T_sat ≈ 179.9 °C.
- Saturation values shown on screen are interpolated from a finely spaced saturation table computed at startup (clustered near the critical point). Small differences in the last printed digit compared with tables are possible.
- The van der Waals option is a qualitative model in reduced units (p/p_c, v/v_c, T/T_c). It reproduces the shape of real behavior but not real numbers. Its saturation curve is computed with the Maxwell equal-area construction.

### Limitations
- **No solid phase.** The surface starts about 1 K above the triple point. Solid, solid-liquid (melting), and solid-vapor (sublimation) regions are not shown, so the p-T diagram shows only the vaporization curve.
- **Logarithmic p and v axes.** These are needed to show the thin liquid region and the huge vapor region on one plot, but they make the surface look different from the stylized linear-axis sketches in most textbooks. It's worth pointing this out to students.
- **Plotted range.** Temperatures run from just above the triple point to 1.4 T_c, and pressures up to 4.5 p_c.
- **Metastable states are not shown.** The van der Waals loops inside the dome (superheated liquid and subcooled vapor) are replaced by the flat equilibrium lines.

### Requirements
Python 3 with numpy, scipy, and matplotlib, plus CoolProp for real-fluid data. Everything can be installed at once with:
```
python -m pip install -r requirements.txt
```
If CoolProp is not available, the program automatically falls back to the van der Waals fluid.

### Feedback and contributions
Suggestions, bug reports, and classroom experiences are welcome through GitHub Issues on this repository. If you adapt the tool for your own course, a note about how you used it helps improve it.
