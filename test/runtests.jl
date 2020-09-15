using Test, ABG

function test_my_styled_printing()
    for i in 1:256
        printstyled(' ', lpad("$i", 3, '0'), ' ', color=i)
        (i % 16 == 0) && println()
    end
    println()

    println("\nAn average messages using color default\n")
    ABG.title("This is a title using bold, color cyan and 154")
    ABG.subtitle("This is a subtitle using color cyan and 154")
    ABG.item("Item 1 using color 74")
    ABG.item("Item 2")
    ABG.message("A description message using color :light_magenta")
    ABG.warning("This is a warning using color 229", "Another warning message")
    ABG.errormsg("One error", "Another error")
    ABG.println()
    ABG.done("Message done, Using color 40.")
    ABG.message("All message tests should be OK")
end

@test 1==1
test_my_styled_printing()
