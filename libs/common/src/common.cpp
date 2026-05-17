#include "common.h"
#include <iostream>
#include <ctime>
#include <iomanip>


namespace orb{

    Logger::Logger(const std::string &file) : filename(file)
    {
        logFile.open(filename, std::ios::app);
        if (!logFile.is_open())
        {
            std::cerr << "Error: could not open log file: " << filename << std::endl;
        }
    }

    Logger::~Logger()
    {
        if (logFile.is_open())
        {
            logFile.close();
        }
    }

    std::string Logger::currentDateTime()
    {
        std::time_t now = std::time(nullptr);
        std::tm *localtm = std::localtime(&now);
        char buffer[80];
        std::strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", localtm);
        return std::string(buffer);
    }

    std::string Logger::levelToString(LogLevel level)
    {
        switch (level)
        {
        case LogLevel::INFO:
            return "INFO";
        case LogLevel::WARNING:
            return "WARNING";
        case LogLevel::ERROR:
            return "ERROR";
        default:
            return "UNKNOWN";
        }
    }

    void Logger::log(LogLevel level, const std::string &message)
    {
        std::string logEntry = "[" + currentDateTime() + "] " + "[" + levelToString(level) + "] " + message;

        // Write to file
        if (logFile.is_open())
        {
            logFile << logEntry << std::endl;
        }

        // ANSI colors
        std::string colorReset = "\033[0m";
        std::string colorCode;

        switch (level)
        {
        case LogLevel::INFO:
            colorCode = "\033[32m";
            break; // Green
        case LogLevel::WARNING:
            colorCode = "\033[33m";
            break; // Yellow
        case LogLevel::ERROR:
            colorCode = "\033[31m";
            break; // Red
        default:
            colorCode = "\033[0m";
            break;
        }

        // Print to terminal with color
        std::cout << colorCode << logEntry << colorReset << std::endl;
    }

    void Logger::info(const std::string &message)
    {
        log(LogLevel::INFO, message);
    }

    void Logger::warning(const std::string &message)
    {
        log(LogLevel::WARNING, message);
    }

    void Logger::error(const std::string &message)
    {
        log(LogLevel::ERROR, message);
    }
}