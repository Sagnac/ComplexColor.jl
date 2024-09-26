import Base: range

function range(i::Interval{:open, :closed}; length::Integer)
    range(; stop = rightendpoint(i), step = width(i)/length, length)
end

range(i::Interval{:open, :closed}, length::Integer) = range(i; length)

function range(i::OpenInterval; length::Integer)
    step = width(i) / (length + 1)
    range(leftendpoint(i) + step; step, length)
end

range(i::OpenInterval, length::Integer) = range(i; length)
