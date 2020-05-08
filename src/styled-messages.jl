# The color plan below is better in terminal dark mode
function title(msg::AbstractString)
    printstyled('\n', msg, '\n', bold = true, color = :cyan)
    printstyled(repeat('=', length(msg)), '\n', color=154)
end

function message(msg::AbstractString)
    printstyled('\n', msg, '\n'; color = :light_magenta)
end

function warning(msg::AbstractString)
    printstyled(msg, '\n'; color=229)
end

function item(it::AbstractString)
    printstyled("\n- $it\n"; color=74)
end

function done(msg::AbstractString = "Done")
    printstyled("\t... $msg\n"; color=40)
end

#for i in 1:256
#    printstyled(' ', lpad("$i", 5, '0'), ' ', color=i)
#    (i % 10 == 0) && println()
#end

function test_my_styled_printing()
    title("This is a title")
    item("Item 1")
    item("Item 2")
    message("A description message")
    warning("This is a warning")
    println("An average messages")
    done()
end
