# Command Line Training Tasks

This folder contains exercises to practice basic command line operations.

## Task 1: Create a File with VS Code

1. Open a new file in VS Code from the terminal:
   ```bash
   code fileA
   ```
2. Write the following content in the file:
   ```
   This is file A
   ```
3. Save the file (Ctrl+S or Cmd+S)
4. Verify the file exists in the current directory:
   ```bash
   ls
   ```

## Task 2: Copy and Edit a File

1. Copy `fileA` to `fileB`:
   ```bash
   cp fileA fileB
   ```
2. Open `fileB` in VS Code:
   ```bash
   code fileB
   ```
3. Add the following line to the file:
   ```
   This is a copy fileB from fileA
   ```
4. Save the file

## Task 3: Create Directories

1. Create two folders using the `mkdir` command:
   ```bash
   mkdir folder_A
   mkdir folder_B
   ```
2. Verify the folders were created:
   ```bash
   ls
   ```

## Task 4: Move Files to Folders

1. Move `fileA` to `folder_A`:
   ```bash
   mv fileA folder_A/
   ```
2. Move `fileB` to `folder_B`:
   ```bash
   mv fileB folder_B/
   ```
3. Verify the files have been moved:
   ```bash
   ls folder_A/
   ls folder_B/
   ```

## Task 5: Search and Display File Contents

1. Use `cat` to display the contents of both files:
   ```bash
   cat folder_A/fileA
   cat folder_B/fileB
   ```
2. Use `grep` to search for the word "copy" in all files:
   ```bash
   grep -r "copy" .
   ```
3. Use `rm -r` to remove folder_B:     
   ```bash
   rm -r folder_B/
   ```

## Bonus Tips

- Use `pwd` to see your current working directory
- Use `cd` to navigate between directories
- Use `rm` to remove files (be careful!)
- Use `rm -r` to remove directories and their contents
- Use `help command` to get more information about a command
