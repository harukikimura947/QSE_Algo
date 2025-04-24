@echo off
REM ユーザーに .c ファイルの名前を聞く（拡張子は入力不要）
set /p source_name=Enter source filename (without .c): 

REM 現在の日時を取得してログファイル名に使用
set datetime=%DATE:~0,4%%DATE:~5,2%%DATE:~8,2%_%TIME:~0,2%%TIME:~3,2%%TIME:~6,2%
set datetime=%datetime: =0%

REM フルパスに変換（見やすいログ用）
set full_c=src\core\%source_name%.c
set full_exe=src\core\%source_name%.exe
set full_log=logs\%source_name%_%datetime%.log

REM ビルドステップ
echo === Compiling %full_c% ...
gcc %full_c% -o %full_exe%
if errorlevel 1 (
    echo ❌ Compilation failed.
    pause
    exit /b
)

REM 実行ステップ
echo === Running %full_exe% ...
cd src\core
%source_name%.exe > ..\..\%full_log% 2>&1

echo === Execution complete. Log saved to %full_log%
pause