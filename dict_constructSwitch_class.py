class Choose_Cont:
    def __init__(self):
        pass

    def PrintBlue():
        return("some people choose Blue!\n") 
    
    def PrintRed():
        return("some people choose Red!\n") 
    
    def PrintOrange():
        return("some people choose Orange!\n")

    def PrintYellow():
        return("some people choose Yellow!\n")

    global switch_flag
    switch_flag = {"Ame" : PrintBlue, "Eur":PrintRed, "Aus":PrintOrange, "Asi":PrintYellow }

    def switch(self, continent):
        return switch_flag.get(continent,  "invalid area.")

a = Choose_Cont()
print(a.switch('Ame')())
#a.switch('Ame')


