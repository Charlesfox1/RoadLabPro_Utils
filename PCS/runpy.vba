Sub runpython()
    Dim Ret_Val
    Dim args As String

    args = "W:\programming\python\other_py\sqrt.py"
    Ret_Val = Shell("C:\Program Files (x86)\python27\python.exe" & " " & args, vbNormalFocus)
    If Ret_Val = 1 Then
       MsgBox "Couldn't run python script!", vbOKOnly
    End If
End Sub
